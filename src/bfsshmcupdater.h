/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef BFSSHMCUPDATER_H
#define BFSSHMCUPDATER_H

#include "bfssbosonicaction.h"
#include "bfsssystemstate.h"
#include "bfssconfigwriter.h"
#include "configparser.h"
#include "RandomNumberGenerator.h"
#include "bfssdiracaction.h"
#include <fstream>

struct HmcParameter{
    Real deltaTauScalar=0.0;
    Real deltaTauGauge=0.0;
    size_t numStepsBosonic=0;
    size_t numStepsFermionic=0;
    bool hotstart=false;
    bool checkhmc=false;
    unsigned long approxCheck=0;
    DiracActionParameter diracParameter={};
};

bool readHmcParameter(Parser::ConfigParser& parser,HmcParameter& parameter);

struct HmcObservables{
    Real oldaction=0.0;
    Real momentumpart=0.0;
    unsigned long accept=0;
    unsigned long reject=0;
    unsigned long counter=0;
    unsigned long alphaconstrcounter=0;
    unsigned long failcounter=0;
    Real deltaH=0;
    Real deltaHSquare=0.0;
    Real expdeltaH=0.0;
    Real expdeltaHSquare=0.0;
    bool fail=false;
    void updateFail(bool in){
      fail=fail||in;
    }
    void reset(){
      accept = 0;
      reject = 0;
      alphaconstrcounter=0;
      deltaH = 0.0;
      deltaHSquare = 0.0;
      expdeltaH = 0.0;
      expdeltaHSquare = 0.0;
      counter = 0;
      fail=false;
    }
    void addDeltaH(const Real dH){
      const Real edH(exp(dH));
      deltaH += dH;
      deltaHSquare += dH * dH;
      expdeltaH += edH;
      expdeltaHSquare += edH * edH;
    }
};



class BfssUpdaterHmc {

  public:
    using Configuration=BfssConfiguration;
    using Mat=BfssConfiguration::MatrixType;
    using Rand=RandomNumberGenerator;
    using SystemState=BfssSystemState;

private:
  HmcParameter parameter_={};
  Rand rnd_={};
	BfssBosonicAction action_={};
	DiracAction* diracaction_=nullptr;
	HmcObservables observables_={};
	Configuration momentum_={};
	Configuration force_={};
	Configuration oldconfig_={};
	std::ofstream hmclog_=std::ofstream("hmc.log",std::ios_base::app);
public:
    
	BfssUpdaterHmc() {

	}
  virtual ~BfssUpdaterHmc() {
    delete diracaction_;
  }

	virtual void reInitialize(SystemState& in) {
		reset();
	}

	virtual void reset() {
	  observables_.reset();
	  hmclog_<<"\n#Reset:1)confignr 2)counter 3)dH 4)oldH 5)newH 6)accept 7)alphaConstraint 8)Fail 9)FermionCGCount"<<std::endl;
	}

	void initialize(BfssSystemState& state,
			Parser::ConfigParser& parser) {
	  force_.resize(state.config().size());
	  momentum_.resize(state.config().size());
	  oldconfig_.resize(state.config().size());
	  reset();
    observables_.updateFail(readHmcParameter(parser,parameter_));
    parameter_.diracParameter.latticesize=state.config().size();
    if(parameter_.diracParameter.rhmc&&!diracaction_){
      diracaction_=new DiracActionRHmc(state,parameter_.diracParameter); // keeps reference to state!
    }
    if(parameter_.diracParameter.hmc&&!diracaction_){
      diracaction_=new DiracActionHmc(state,parameter_.diracParameter); // keeps reference to state!
    }
	}

	void initializeFields(SystemState& state) {
	  if(parameter_.hotstart)
     setHotStart(state.config());
	  else
	    setColdStart(state.config());
	}

	void runChecks(SystemState& state) {
	  if(parameter_.checkhmc){
	   checkForce(state);
	   checkIntegrator(state);
	  }
	}



	void update(SystemState& in) {
		oldconfig_ = in.config();
		bool fail=false;

    initializeMomenta();
    if(diracaction_){
      if(parameter_.approxCheck>0 && observables_.counter%parameter_.approxCheck==0){
        diracaction_->setCheck(true);
      }
      diracaction_->updateNoiseVector(rnd_);
      if(parameter_.approxCheck>0 && observables_.counter%parameter_.approxCheck==0){
        LOG_STREAM<<"Dirac approx check "<<diracaction_->getCheckValue()<<std::endl;
        diracaction_->setCheck(false);
      }
    }
    observables_.oldaction=completeAction(in);
		integrateLeapFrog(in);
		const Real newa(completeAction(in));
		const Real deltaH(newa - observables_.oldaction);
		hmclog_<<in.configNum()<<" "<<observables_.counter<< " " << deltaH << " " << observables_.oldaction << " " << newa;
		if (newa != newa) {
			LOG_MESSAGE("Error received NAN values for new action");
			action_.logActionParts(LOG_STREAM);
			LOG_STREAM<<std::endl;
			fail=true;
		}
		fail=fail||checkFailedDiracAction();
		bool accept=false;
		bool constr=false;
		if(!fail){
    const Real r(Real(1.0) - rand()); //[0,1) to (0,1]
		constr=constraintViolationAlpha(in.config());
		if(constr){
		  observables_.alphaconstrcounter++;
		}
		accept=(log(r) <= -deltaH);
#ifdef ACCEPT_FIRST
		if(in.configNum()<10){
		  accept=true;
		}
#endif
		if (accept&&!constr) {
			 observables_.accept++;
			 observables_.oldaction=newa;
	     projectU1Part(in.config());
	     notifyFieldChange(); // since projection is done.
		} else {
			in.config() = oldconfig_;
			notifyFieldChange();
			observables_.reject++;
		}
    if(!constr){
      observables_.addDeltaH(deltaH);
    }
		}else{ // failed
      observables_.failcounter++;
      in.config() = oldconfig_;
      notifyFieldChange();
		}
    hmclog_ <<" "<<accept<<" "<<constr<<" "<<fail<<" "<<diracActionIterationCount();
		hmclog_<<std::endl;
		observables_.counter++;
	}


  void checkIntegrator(SystemState& in){
    std::cout << "Checking integrator" << std::endl;
    initializeMomenta();

    if(diracaction_)
      diracaction_->updateNoiseVector(rnd_);

    const auto oldnumstepsb=parameter_.numStepsBosonic;
    const auto oldnumstepsf=parameter_.numStepsFermionic;
    const auto olddtalpha=parameter_.deltaTauGauge;
    const auto olddts=parameter_.deltaTauScalar;

    const Real refaction=completeAction(in);
    const Real refactionferm=fermionAction(in);
    std::ofstream ofs("integrator_check.log");
    ofs<<"#dTauG dTauS dH reversDH Hfinal reversDHFermion"<<std::endl;
    for(int i(0);i<10;i++) {
      integrateLeapFrog(in);
      const Real actionInt=completeAction(in);
      parameter_.deltaTauScalar*=-1;
      parameter_.deltaTauGauge*=-1;
      integrateLeapFrog(in);
      parameter_.deltaTauScalar*=-1;
      parameter_.deltaTauGauge*=-1;
      const Real actionIntBW=completeAction(in);
      const Real actionIntBWFerm=fermionAction(in);
      ofs<<parameter_.deltaTauScalar<<" "<<parameter_.deltaTauGauge<<" "<<actionInt-refaction<<" "<<actionIntBW-refaction<<" "<<actionInt<<" "<<actionIntBWFerm-refactionferm<<std::endl;
      parameter_.numStepsBosonic*=1.5;
      parameter_.deltaTauScalar=olddts*Real(oldnumstepsb)/Real(parameter_.numStepsBosonic);
      parameter_.deltaTauGauge=olddtalpha*Real(oldnumstepsb)/Real(parameter_.numStepsBosonic);
    }

    parameter_.numStepsBosonic=oldnumstepsb;
    parameter_.numStepsFermionic=oldnumstepsf;
    parameter_.deltaTauGauge=olddtalpha;
    parameter_.deltaTauScalar=olddts;
  }


	void checkForce(SystemState& in) {
		LOG_MESSAGE("Checking gauge and complete force");
		std::ofstream ofs("force_check.log");
		ofs<<"#delta tau : force - force_approx : force_approx : force : constraint violation\n";

		initializeMomenta();
		force_.setToZero();
		action_.addForce(in,force_,1.0);
    if(diracaction_){
      diracaction_->updateNoiseVector(rnd_);
      diracaction_->addForce(force_,1.0);
    }
		const Real refvalue=2.0*(sumProductTraceX(momentum_,force_)
		    +sumProductAlpha(momentum_,force_));
		const Real dtaurange(parameter_.deltaTauScalar/Real(10));
    Real deltatauref = dtaurange;
		while (deltatauref > dtaurange/100.0) {
			fieldStep(in,-0.5*deltatauref,-0.5*deltatauref);
			const Real newaction1(completeAction(in));
			bool constr=constraintViolationAlpha(in.config());
      fieldStep(in,deltatauref,deltatauref);
      const Real newaction2(completeAction(in));
      fieldStep(in,-0.5*deltatauref,-0.5*deltatauref);
			const Real forceappox((newaction2 - newaction1) / deltatauref);
			ofs << deltatauref << " " << fabs(refvalue - forceappox) << " "
					<< forceappox << " " << refvalue<<" "<<constr<< std::endl;
			deltatauref -= dtaurange/10.0;
		}
	}

  void setHotStart(Configuration& in) {
    setGaussianRandomX(in,rnd_);
    setUniformRandomAlpha(in,rnd_);
    reduceToFitConstraint(in);
    projectU1Part(in);
    notifyFieldChange();
  }

  void setColdStart(Configuration& in) {
    in.setPhaseToIdentity();
    in.setToZero();
//    setFortranCodeAlpha(in);
    setUniformRandomAlpha(in,rnd_);
    reduceToFitConstraint(in);
    projectU1Part(in);
    notifyFieldChange();
  }

  void projectU1Part(Configuration& in)const{
    traceOutAlpha(in);
    traceOutX(in);
  }

  void printData(std::ostream& os,const SystemState& state)const{
    os<<"Finished "<<observables_.counter<<"  updates\n";
    os<<"Acceptance rate: "<<double(observables_.accept)/double(observables_.counter)<<"\n";
    os<<"Alpha constraint violations "<<observables_.alphaconstrcounter<<"\n";
    os<<"Failed HMC updates "<<observables_.failcounter<<"\n";
    const double r(1.0/double(observables_.counter));
    os<<"Average deltaH "<<observables_.deltaH*r<<" "<<std::sqrt(r*(observables_.deltaHSquare*r-observables_.deltaH*observables_.deltaH*r*r))<<std::endl;
    os<<"Average exp(deltaH) "<<observables_.expdeltaH*r<<" "<<std::sqrt(r*(observables_.expdeltaHSquare*r-observables_.expdeltaH*observables_.expdeltaH*r*r))<<std::endl;
  }

private:

  bool checkFailedDiracAction(){
    if(diracaction_){
      if(diracaction_->fail()){
        LOG_MESSAGE("Error: Failed inverter in Dirac action")
        return true;
      }
    }
    return false;
  }

  unsigned long diracActionIterationCount(){
    if(diracaction_){
        return diracaction_->iterationCount();
    }
    return 0;
  }

  void notifyFieldChange(){
    if(diracaction_)
      diracaction_->setGaugeFieldDirty();
  }

	Real rand() {
		return rnd_();
	}

  void initializeMomenta() {
    setGaussianRandomX(momentum_,rnd_);
    setGaussianRandomAlpha(momentum_,rnd_);
  }



	void integrateLeapFrog2(SystemState& state) {
		if (parameter_.numStepsBosonic == 0)
			return;
		const Real dtx=parameter_.deltaTauScalar;
		const Real dtalpha=parameter_.deltaTauGauge;
		momentumStep(state, 0.5 * dtx,0.5*dtalpha);
		for (size_t i(parameter_.numStepsBosonic - 1); i--;) {
		  fieldStep(state,dtx,dtalpha);
		  momentumStep(state,dtx,dtalpha);
		}
		fieldStep(state,dtx,dtalpha);
		momentumStep(state, 0.5 * dtx,0.5*dtalpha);
	}

  void integrateLeapFrog(SystemState& state) {
    if (parameter_.numStepsBosonic == 0)
      return;
    const Real dtx=parameter_.deltaTauScalar;
    const Real dtalpha=parameter_.deltaTauGauge;
    fieldStep(state, 0.5 * dtx,0.5*dtalpha);
    for (size_t i(parameter_.numStepsBosonic - 1); i--;) {
      momentumStep(state,dtx,dtalpha);
      fieldStep(state,dtx,dtalpha);
    }
    momentumStep(state,dtx,dtalpha);
    fieldStep(state, 0.5 * dtx,0.5*dtalpha);
  }

	Real momentumPart() {
		return sumTraceX2(momentum_)+sumAlpha2(momentum_);
	}

	Real bosonicAction(const SystemState& in) {
		action_.update(in);
		return action_.action();
	}

  Real fermionAction(const SystemState& in) {
    if(diracaction_){
      diracaction_->update();
      return diracaction_->action();
    }
    return 0;
  }

	Real completeAction(const SystemState& state) {
		const Real bact(bosonicAction(state));
		const Real mpart(momentumPart());
		const Real faction(fermionAction(state));
		return mpart + bact+faction;
	}
	
	void momentumStep(const SystemState& state, const Real dtx,const Real dtalpha) {
	  momentumStepComplete(state,dtx,dtalpha);
	}

	void momentumStepBosonic(const SystemState& state, const Real dtx,const Real dtalpha) {
	    force_.setToZero();
	    action_.addForce(state,force_,1.0);
	    addRescaled(momentum_,force_,-dtx,-dtalpha);
	}

  void momentumStepFermonic(const SystemState& state, const Real dtx,const Real dtalpha) {
      force_.setToZero();
      if(diracaction_)
        diracaction_->addForce(force_,1.0);
      addRescaled(momentum_,force_,-dtx,-dtalpha);
  }

  void momentumStepComplete(const SystemState& state, const Real dtx,const Real dtalpha) {
      force_.setToZero();
      action_.addForce(state,force_,1.0);
      if(diracaction_)
        diracaction_->addForce(force_,1.0);
      addRescaled(momentum_,force_,-dtx,-dtalpha);
  }


	void fieldStep(SystemState& in, const Real dtx,const Real dtalpha) {
	  if(in.parameters().backroundgauge){
	    addRescaled(in.config(),momentum_,dtx,0.0);
	  }else{
     addRescaled(in.config(),momentum_,dtx,dtalpha);
	  }
	  notifyFieldChange();
	}


};


inline bool readHmcParameter(Parser::ConfigParser& parser,HmcParameter& parameter){
  Real dtau=0.0;
  bool fail=false;
  try{
    Real tautotal=parser.getdouble("HMC","TauTotal");
    unsigned int nsteps=parser.getdouble("HMC","NumSteps");
    dtau=tautotal/Real(nsteps);
  }catch(...){
  }
  try{
    dtau=parser.getdouble("HMC","DTau");
  }catch(...){
  }
  parameter.deltaTauGauge=dtau;
  parameter.deltaTauScalar=dtau;
  try{
    parameter.deltaTauScalar=parser.getdouble("HMC","DTauScalar");
  }catch(...){
  }
  try{
    parameter.deltaTauGauge=parser.getdouble("HMC","DTauGauge");
  }catch(...){
  }
  if(parameter.deltaTauGauge==0.0 || parameter.deltaTauScalar==0.0){
    LOG_STREAM<<"Error: parameter DTau or DTauScalar and DTauGauge not set in Section HMC!"<<std::endl;
    fail=true;
  }
  try{
    parameter.numStepsBosonic=parser.getdouble("HMC","NumSteps");
  }catch(...){
    LOG_STREAM<<"Error: parameter NumSteps not set in Section HMC!"<<std::endl;
    fail=true;
  }
  try{
    parameter.numStepsFermionic=parser.getdouble("HMC","NumStepsFermion");
  }catch(...){
    LOG_STREAM<<"Parameter NumStepsFermion not set in Section HMC using NumSteps"<<std::endl;
    parameter.numStepsFermionic=parameter.numStepsBosonic;
  }
  parameter.hotstart=parser.get_value_noerr("HMC","HotStart",true);
  parameter.checkhmc=parser.get_value_noerr("HMC","CheckHMC",false);
  parameter.approxCheck=parser.get_value_noerr("HMC","RHMCCheck",0);
  parameter.diracParameter.hmc=parser.get_value_noerr("HMC","FHMC",false);
  parameter.diracParameter.rhmc=parser.get_value_noerr("HMC","FRHMC",false);
  parameter.diracParameter.residuum=parser.get_value_noerr("HMC","Resid",Real(1e-4));
  parameter.diracParameter.maxIterations=parser.get_value_noerr("HMC","MaxIt",1000);
  try{
   parameter.diracParameter.filenameActionRat=parser.get("HMC","ActionRat");
   parameter.diracParameter.filenameForceRat=parser.get("HMC","ForceRat");
  }catch(...){

  }
  parameter.diracParameter.rescaleapprox=parser.get_value_noerr("HMC","RatApproxLambdaMax",parameter.diracParameter.rescaleapprox);
  LOG_STREAM<<"Parameter read:DTauScalar:"<<parameter.deltaTauScalar<<":DTauGauge:"<<parameter.deltaTauGauge
      <<":NumSteps:"<<parameter.numStepsBosonic<<":NumStepsFermion:"<<parameter.numStepsFermionic
      <<":HotStart:"<<parameter.hotstart<<":CheckHMC:"<<parameter.checkhmc<<":FHMC:"<<parameter.diracParameter.hmc<<":FRHMC:"<<parameter.diracParameter.rhmc<<std::endl;
  LOG_STREAM<<"Fermion action:inversion residuum(Resid):"<<parameter.diracParameter.residuum<<":MaxIt:"<<parameter.diracParameter.maxIterations<<std::endl;
  if(parameter.diracParameter.rhmc){
    LOG_STREAM<<"Fermion action rhmc:ActionRat:"<<parameter.diracParameter.filenameActionRat<<":ForceRat:"
        <<parameter.diracParameter.filenameForceRat<<":RHMCCheck:"<<parameter.approxCheck<<":RatApproxLambdaMax:"<<parameter.diracParameter.rescaleapprox<<std::endl;
  }
  return !fail;
}

#endif
