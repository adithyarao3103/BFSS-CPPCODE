
/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 *
 *  Adjusted implementation from older SYMQM and SYM project.
 */


#ifndef SYMQMDIRACFORCE_H
#define SYMQMDIRACFORCE_H
#include "bfssconfig.h"
#include "diracoperator.h"
#include "cgmsolver.h"
#include "bfsssystemstate.h"
#include "configparser.h"

//#define CHECK_ACTION_MM

struct DiracActionParameter{
    bool hmc=false;
    bool rhmc=true;
    Real residuum=1e-7;
    unsigned int maxIterations=1000;
    size_t latticesize=0;
    std::string filenameActionRat={};
    std::string filenameForceRat={};
    Real rescaleapprox=1.0;
};

class DiracAction {
public:
    using State=BfssSystemState;
    using Config=BfssConfiguration;
    using Vector=BfssDiracVector;
private:
    bool fail_=false;
    bool check_=false;
    Real checkValue_=0.0;
    unsigned long iterationCounter_=0;
public:
	virtual ~DiracAction() {
	}
	virtual void updateNoiseVector(RandomNumberGenerator& rand)=0;
	virtual void update()=0;
  virtual Real action() const=0;
	virtual void addForce(Config &force, const Real prefact)=0;

	virtual void setGaugeFieldDirty()=0;


    bool fail() const {
      return fail_;
    }

    unsigned long iterationCount()const{
      return iterationCounter_;
    }

    void setCheck(const bool check=false){
      check_=check;
    }

    Real getCheckValue()const{
      return checkValue_;
    }

protected:
    bool isCheckMode()const{
      return check_;
    }
    void setCheckValue(const Real v){
      checkValue_=v;
    }
    void setFail(const bool fail = true) {
      fail_ = fail;
    }
    void updateSolverFail(const bool fail){
      fail_=fail_||fail;
    }
    void setIterationCount(const unsigned long count=0){
      iterationCounter_=count;
    }
    void addIterationCount(const unsigned long count){
      iterationCounter_+=count;
    }
};


/**
 * HMC Dirac Action, not needed since applies only for two flavour case, but it helps as a Test.
 */
class DiracActionHmc: public DiracAction {
  public:
      using State=DiracAction::State;
      using Config=DiracAction::Config;
      using Vector=DiracAction::Vector;
      using DiracOperator=BfssDiracOperator;
private:
	Real action_={};
	Vector pf_={};
	Vector chi_={};
	Vector X_={};
	Vector Y_={};
	CgSqSolver<DiracOperator,DiracOperator> solver_;
	bool actiondirty_=true;
  bool phidirty_=true;
	bool Xdirty_=true;
	bool Ydirty_=true;
  DiracActionParameter param_={};
  void setDirty() {
    phidirty_ = Xdirty_ = Ydirty_ = true;
  }
public:

	DiracActionHmc(State& in, const DiracActionParameter& p) :
			action_(0.0), pf_(p.latticesize), chi_(p.latticesize), X_(p.latticesize), Y_(p.latticesize),solver_(in.getDiracOperator(),in.getDiracOperatorDagger()),param_(p){
	  solver_.setResidGoal(p.residuum);
	  solver_.setMaxIterations(p.maxIterations);
	}

	virtual ~DiracActionHmc() {
	}

	void updateNoiseVector(RandomNumberGenerator& rand)override final {
	  setGaussianRandomVector(chi_,rand);
		setDirty();
		action_ = normVectorSq(chi_);
		actiondirty_=false;
		DiracAction::setFail(false);
		DiracAction::setIterationCount();
	}

	void setGaugeFieldDirty()override final {
    actiondirty_=true;
		Xdirty_ = Ydirty_ = true;
	}

	void updatePhi() {
		if (!phidirty_)
			return;
		solver_.opDag().multiply(chi_, pf_); // \phi=D^\dag\chi
		phidirty_ = false;
	}

	void updateX() {
		if (!Xdirty_)
			return;
		solver_.multiply(pf_, X_); //X_=(D^dag D)^-1 pf_
		DiracAction::addIterationCount(solver_.iterations());
    DiracAction::updateSolverFail(solver_.fail());
		LOG_STREAM<<"DiracHMC CG count "<<solver_.iterations()<<std::endl;
		Xdirty_ = false;
	}

	void updateY() {
		if (!Ydirty_)
			return;
		solver_.op().multiply(X_, Y_); //Y_=D X_=(D^\dag)^-1 pf_
		Ydirty_ = false;
	}

	void update() override final{
	  if(!actiondirty_)
	    return;
		updatePhi();
		updateX();
		action_ = real(scalarProd(pf_, X_));
		actiondirty_=false;
	}

	Real action() const override final{
		return action_;
	}

	void addForce(Config &force, const Real prefact)override final{
		updatePhi();
		updateX();
		updateY();
		solver_.op().addForce(force, X_, Y_,prefact);
	}
};



class RationalApproximation {
public:
	typedef std::vector<Real> CoeffVector;
private:
	CoeffVector shift_;
	CoeffVector prefact_;
	int numerator_=1;
	int denominator_=1;
public:

	void setFraction(const int numerator, const int denominator){
	  numerator_=numerator;
	  denominator_=denominator;
	}

	RationalApproximation() :
			shift_(), prefact_() {

	}
	RationalApproximation(const CoeffVector& shift, const CoeffVector& prefact) :
			shift_(shift), prefact_(prefact) {
		if (prefact_.size() == 0) {
			return;
		}
		if (prefact_.size() != shift_.size() + 1) {
			shift_.resize(prefact_.size() - 1);
		}
	}

	/*
	 * x^p=r^{-p} (r x)^p
	 * r^{-p}\alpha_0+\sum_i \frac{r^{-1-p}\alpha_i}{x-\beta_i/r}
	 * r ~ \lambda_max^{-1}
	 */
	void rescale(const Real fact){
	  const Real power(Real(numerator_)/Real(denominator_));
    prefact_[0]*=std::pow(fact,-power);
    for (size_t i(1); i < prefact_.size(); i++) {
      prefact_[i]*=std::pow(fact,-1.0-power);
      shift_[i-1]/=fact;
    }

	}

	RationalApproximation(const CoeffVector& completerat) :
			shift_(), prefact_() {
		resetPrefactShift(completerat);
	}

	void resetPrefactShift(const CoeffVector& completerat) {
		if (completerat.size() == 0)
			return;
		if (completerat.size() % 2 == 0)
			return;
		size_t npref((completerat.size() - 1) / 2);
		shift_.resize(0);
		prefact_.resize(0);
		for (size_t i(0); i < completerat.size(); i++) {
			if (i > npref) {
				shift_.push_back(completerat[i]);
			} else {
				prefact_.push_back(completerat[i]);
			}
		}
	}

	Real operator()(const Real in) const {
		if (prefact_.size() == 0) {
			return (0.0);
		}
		Real ret(prefact_[0]);
		for (size_t i(1); i < prefact_.size(); i++) {
			ret += prefact_[i] / (in + shift_[i - 1]);
		}
		return (ret);
	}

	const CoeffVector& prefact() const {
		return prefact_;
	}

	void setPrefact(const CoeffVector& prefact) {
		prefact_ = prefact;
	}

	const CoeffVector& shift() const {
		return shift_;
	}

	void setShift(const CoeffVector& shift) {
		shift_ = shift;
	}
	size_t shiftSize() const {
		return (shift_.size());
	}

	void print(std::ostream& os) {
		if (prefact_.size() == 0) {
			os << "#no rational approx";
			return;
		}
		os << prefact_[0] << "\n";
		for (size_t i(1); i < prefact_.size(); i++) {
			os << prefact_[i] << "," << shift_[i - 1] << "\n";
		}
		os<<"approximation of x^"<<numerator_<<"/"<<denominator_<<"\n";
	}

};

void readRationalApproximation(RationalApproximation& approx, const std::string& filename,
    const int ind) {
  const std::string alname("alpha");
  const std::string bename("beta");
  std::ifstream file(filename.c_str());
  std::vector<Real> prefact;
  std::vector<Real> shifts;
  int counter(0);
  bool poly(false);
  while (file.good()) {
    std::string line;
    bool firstitem(false);
    std::getline(file, line);
    if (line.find(alname) != std::string::npos) {
      if (poly == false) {
        counter++;
        poly = true;
        firstitem = true;
      }
    } else {
      poly = false;
    }
    if (counter == ind + 1) {
      if (poly) {
        std::string num(
            line.substr(line.find("[") + 1, line.find("]")));
        std::string value(line.substr(line.find("=") + 1));
        std::istringstream numstream(num);
        std::istringstream valstream(value);
        int ind(0);
        Real val(0);
        numstream >> ind;
        valstream >> val;
        prefact.push_back(val);
      }
      if (poly && !firstitem) {
        std::size_t betapos(line.find("beta"));
        std::string num(
            line.substr(line.find("[", betapos) + 1,
                line.find("]", betapos)));
        std::string value(line.substr(line.find("=", betapos) + 1));
        std::istringstream numstream(num);
        std::istringstream valstream(value);
        int ind(0);
        Real val(0);
        numstream >> ind;
        valstream >> val;
        shifts.push_back(val);
      }
    }
  }
  approx.setPrefact(prefact);
  approx.setShift(shifts);

}


class  DiracActionRHmc: public DiracAction {
public:
    using State=BfssSystemState;
    using Config=BfssConfiguration;
    using Vector=BfssDiracVector;
	using VectorVector=std::vector<Vector>;
  using DiracOperator=BfssDiracOperator;
private:
	Real action_=0.0;
	Vector chi_={}; //< random noise
	Vector pf_={}; //< pseudofermion
	VectorVector tmpphi_={}; //< for calculation of pseudofermion
	RationalApproximation approxphi_={}; // approximation of x^Nf/2
  RationalApproximation approxforce_={};// approximation of x^-Nf/4
	VectorVector X_={};
	VectorVector Y_={};
	CgSqMSolver<DiracOperator,DiracOperator> solver_;
  bool actiondirty_=true;
  bool phidirty_=true;
  bool Xdirty_=true;
  bool Ydirty_=true;
  DiracActionParameter param_={};
  void setDirty() {
    phidirty_ = Xdirty_ = Ydirty_ = true;
  }

  void adjustVectorVector(VectorVector& in, const size_t size)const{
    for(auto& c:in){
      c.resize(size);
    }
  }


public:
  DiracActionRHmc(State& in, const DiracActionParameter& p) :
    action_(0.0), pf_(p.latticesize), chi_(p.latticesize),solver_(in.getDiracOperator(),in.getDiracOperatorDagger()),param_(p){
    solver_.setResidGoal(p.residuum);
    solver_.setMaxIterations(p.maxIterations);
    setRationalApprox(p);
 }

	void setRationalApprox(const DiracActionParameter& p) {
    readRationalApproximation(approxphi_,p.filenameActionRat,0);
    approxphi_.setFraction(1,8);
    readRationalApproximation(approxforce_,p.filenameForceRat,1);
    approxforce_.setFraction(-1,4);
    LOG_STREAM << "Using the following rational approx Action:\n";
		approxphi_.print(LOG_STREAM);
    if(approxphi_.shiftSize()==0||approxphi_.shift().size()+1!=approxphi_.prefact().size()){
     LOG_STREAM<<"Error with phi approx"<<std::endl;
     throw;
    }
    LOG_STREAM << "Using the following rational approx Force:\n";
		approxforce_.print(LOG_STREAM);
    if(approxforce_.shiftSize()==0||approxforce_.shift().size()+1!=approxforce_.prefact().size()){
     LOG_STREAM<<"Error with force approx"<<std::endl;
     throw;
    }
    LOG_STREAM << "Using the following rescaled rational approx Action:\n";
    approxphi_.rescale(1.0/param_.rescaleapprox);
    approxphi_.print(LOG_STREAM);
    LOG_STREAM << "Using the following rescaled rational approx Force:\n";
    approxforce_.rescale(1.0/param_.rescaleapprox);
    approxforce_.print(LOG_STREAM);
		tmpphi_.resize(approxphi_.shiftSize());
		X_.resize(approxforce_.shiftSize());
		Y_.resize(approxforce_.shiftSize());
    adjustVectorVector(tmpphi_,pf_.size());
    adjustVectorVector(X_,pf_.size());
    adjustVectorVector(Y_,pf_.size());
		testParameters();
    solver_.setShiftVector(approxforce_.shift());
	}

	void testParameters() {
	  LOG_STREAM << "Testing rational parameters" << std::endl;
		std::ofstream of("rationalparameters_check.log");
		double fakt(1e-9);
		for (int j(10); j--;) {
			double x = fakt;
			for (int i(0); i < 1000; i++) {
				const Real x18(approxphi_(x));
				const Real xm14(approxforce_(x));
				of << x << " " << x18 << " " << xm14 << " "
						<<x18-std::pow(x,1.0/8.0)<<" "
						<<xm14-std::pow(x,-1.0/4.0)<<" "
						<< x18 * x18 - 1.0 / xm14 << "\n";
				x += fakt;
			}
			fakt *= 10.0;
		}
	}

  void updateNoiseVector(RandomNumberGenerator& rand)override final {
    setGaussianRandomVector(chi_,rand);
    setDirty();
    action_ = normVectorSq(chi_);
    actiondirty_=false;
    DiracAction::setFail(false);
    DiracAction::setIterationCount();
    if(DiracAction::isCheckMode()){
      const Real act(action_);
      actiondirty_=true;
      update();
      DiracAction::setCheckValue(fabs(action_-act)/(action_+act));
      action_=act;
    }
  }


  void setGaugeFieldDirty()override final {
    actiondirty_=true;
    Xdirty_ = Ydirty_ = true;
  }

  void sumUpRational(const RationalApproximation& rat, const VectorVector& vv,
      const Vector& first, Vector& out) const {

#pragma omp parallel for
    for (size_t site = 0; site < out.size(); site++) {
      for (size_t c(out.elementnum); c--;) {
        out(site,c) = rat.prefact()[0] * first(site,c);
      }
    }

    for (size_t k(0); k < rat.shiftSize(); k++) {
#pragma omp parallel for
      for (size_t site = 0; site < out.size(); site++) {
        for (size_t c(out.elementnum); c--;) {
          out(site,c) += rat.prefact()[k + 1] * vv[k](site,c);
        }
      }
    }
  }

	void updatePhi() {
		if (!phidirty_)
			return;
		solver_.setShiftVector(approxphi_.shift());
		solver_.multiply(chi_, tmpphi_);
		solver_.setShiftVector(approxforce_.shift());
		sumUpRational(approxphi_, tmpphi_, chi_, pf_);
    DiracAction::addIterationCount(solver_.iterations());
    DiracAction::updateSolverFail(solver_.fail());
    LOG_STREAM<<"DiracRHMC action CG count "<<solver_.iterations()<<std::endl;
		phidirty_ = false;
	}



	void updateX() {
		if (!Xdirty_)
			return;
		solver_.multiply(pf_, X_);
    DiracAction::addIterationCount(solver_.iterations());
    DiracAction::updateSolverFail(solver_.fail());
    LOG_STREAM<<"DiracRHMC force CG count "<<solver_.iterations()<<std::endl;
		Xdirty_ = false;
	}


	void updateY() {
		if (!Ydirty_)
			return;
		for (size_t i(X_.size()); i--;) {
		  solver_.op().multiply(X_[i], Y_[i]);
		}
		Ydirty_ = false;
	}

	void update() override final{
    if(!actiondirty_)
      return;
    updatePhi();
    updateX();
    sumUpRational(approxforce_, X_, pf_, tmpphi_[0]);
    action_ = real(scalarProd(pf_, tmpphi_[0]));
    actiondirty_=false;
	}

	Real action() const override final{
		return (action_);
	}

  void addForce(Config &force, const Real prefact)override final{
    updatePhi();
    updateX();
    updateY();
    for (size_t k(0); k < X_.size(); k++) {
     solver_.op().addForce(force, X_[k], Y_[k],prefact*approxforce_.prefact()[k + 1]);
    }
  }

};


#endif
