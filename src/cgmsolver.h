/*
 * cgmsolver.h
 *
 *  Created on: 30.03.2023
 *      Author: Georg Bergner
 *
 *  This contains basic solvers, adapted from code developed for the supersymmetric Yang-Mills project.
 */

#ifndef SRC_BFSS_CGMSOLVER_H_
#define SRC_BFSS_CGMSOLVER_H_

//#define DEBUG_INVERTER
#include "bfssconfig.h"
#include <vector>
#include <iostream>

template<class OP, class OPDAG>
class CgnrSolver {
  public:
    using Vector=typename OP::Vector;
  protected:
    OP op_;
    OPDAG opdag_;
    Real residgoal_ = 1e-4;
    Real residsq_ = 0.0;
    unsigned long maxit_ = 1000;
    bool fail_ = false;
    unsigned int project_ = 100;
    unsigned long iterations_ = 0;
  public:
    CgnrSolver(const OP &op, const OPDAG &opdag) :
        op_(op), opdag_(opdag) {
    }

    OP& op(){
      return op_;
    }

    OPDAG& opDag(){
      return opdag_;
    }

    void setResidGoal(const Real res) {
      residgoal_ = res;
    }
    void setMaxIterations(const unsigned long mit) {
      maxit_ = mit;
    }
    bool fail() const {
      return fail_;
    }
    Real residuum() const {
      return sqrt(residsq_);
    }

    unsigned long iterations()const{
      return iterations_;
    }

    void multiply(const Vector &in, Vector &out) {
      const Real normInSq(normVectorSq(in));
      fail_ = false;
      iterations_ = 0;
      Vector rv(in.size());
      Vector dv(in.size());
      Vector tmp(in.size());
      op_.multiply(out, rv);
      rescaleAdd(rv, -1.0, in); //rv=in-rv;
      residsq_ = normVectorSq(rv);
      dv.setToZero();
      opdag_.multiply(rv, dv);
      Real delta = normVectorSq(dv);
      std::cout << "Checking trace 0 " << sumAbsTraceVector(dv) << " " << sumAbsTraceVector(out) << " " << sumAbsTraceVector(rv) << std::endl;
      for (size_t ite = 0; ite < maxit_; ++ite) {
#ifdef DEBUG_INVERTER
        std::cout << "CGNR:" << ite << " " << sqrt(residsq_) <<" "<<normVectorSq(dv)<<" "<<normVectorSq(rv)<< " " << residgoal_ <<" "<<normInSq<< std::endl;
#endif

        if ((residsq_ < residgoal_ * residgoal_ * normInSq || residsq_ != residsq_)) {
          if (residsq_ != residsq_) {
            std::cout << "CGNRERR:" << residsq_ << std::endl;
            out = in;
          }
          op_.multiply(out, rv);
          rescaleAdd(rv, -1.0, in);
          const Real realresidsq(normVectorSq(rv));
#ifdef DEBUG_INVERTER
          std::cout << "CGNR:" << ite << " " << sqrt(residsq_) << " " << sqrt(realresidsq) << " " << residgoal_<<" "<<normInSq << std::endl;
#endif
          residsq_ = realresidsq;
          if (residsq_ < residgoal_ * residgoal_ * normInSq) {
#ifdef DEBUG_INVERTER
            std::cout << "CGNR:Converged:" << ite << " " << sqrt(residsq_) << " " << sqrt(realresidsq) << " " << residgoal_ <<" "<<normInSq<< std::endl;
#endif
            iterations_ = ite;
            return;
          }
          opdag_.multiply(rv, dv);
          delta = normVectorSq(dv);
        }

        op_.multiply(dv, tmp);
        const Real alpha = delta / normVectorSq(tmp);
        addRescaled(out, dv, alpha);
        addRescaled(rv, tmp, -alpha);
        residsq_ = normVectorSq(rv);
        opdag_.multiply(rv, tmp);
        const Real deltanew = normVectorSq(tmp);
        const Real beta = deltanew / delta;
        delta = deltanew;
        rescaleAdd(dv, beta, tmp);
        if (project_ != 0 && ite % project_ == 0) {
          projectTracelessVector(out);
          projectTracelessVector(dv);
          projectTracelessVector(rv);
        }
      }
      std::cout << "CGNR did not converge " << sqrt(residsq_) << " " << residgoal_ << std::endl;
      iterations_ = maxit_;
      fail_ = true;
    }
};
template<class OP, class OPDAG>
class CgSqSolver: public CgnrSolver<OP, OPDAG> {
  public:
    using TBase=CgnrSolver<OP,OPDAG>;
    using Vector=typename TBase::Vector;

    CgSqSolver(const OP &op, const OPDAG &opdag) :
        TBase(op, opdag) {
    }

    void multiply(const Vector &in, Vector &out) {
      const Real normInSq(normVectorSq(in));
      TBase::iterations_ = 0;
      TBase::fail_ = false;
      Vector rv(in.size());
      Vector dv(in.size());
      Vector tmp(in.size());
      Vector tmp2(in.size());
      TBase::op_.multiply(out, tmp);
      TBase::opdag_.multiply(tmp, rv);
      rescaleAdd(rv, -1.0, in); //rv=in-rv;
      dv = rv;
      TBase::residsq_ = normVectorSq(rv);
      for (size_t ite = 0; ite < TBase::maxit_; ++ite) {
#ifdef DEBUG_INVERTER
        std::cout << "CGSQ:" << ite << " " << sqrt(TBase::residsq_) <<" "<<normVectorSq(dv)<<" "<<normVectorSq(rv)<< " " << TBase::residgoal_ <<" "<<normInSq<< std::endl;
#endif

        if ((TBase::residsq_ < TBase::residgoal_ * TBase::residgoal_ * normInSq || TBase::residsq_ != TBase::residsq_)) {
          if (TBase::residsq_ != TBase::residsq_) {
            std::cout << "CGSQERR:" << TBase::residsq_ << std::endl;
            out = in; // try reset.
          }
          TBase::op_.multiply(out, tmp);
          TBase::opdag_.multiply(tmp, rv);
          rescaleAdd(rv, -1.0, in);
          const Real realresidsq(normVectorSq(rv));
#ifdef DEBUG_INVERTER
          std::cout << "CGNR:" << ite << " " << sqrt(TBase::residsq_) << " " << sqrt(realresidsq) << " " << TBase::residgoal_<<" "<<normInSq << std::endl;
#endif
          TBase::residsq_ = realresidsq;
          if (TBase::residsq_ < TBase::residgoal_ * TBase::residgoal_ * normInSq) {
#ifdef DEBUG_INVERTER
            std::cout << "CGSQ:Converged:" << ite << " " << sqrt(TBase::residsq_) << " " << sqrt(realresidsq) << " " << TBase::residgoal_ <<" "<<normInSq<< std::endl;
#endif
            TBase::iterations_ = ite;
            return;
          }
          dv = rv;
        }

        TBase::op_.multiply(dv, tmp);
        TBase::opdag_.multiply(tmp, tmp2);
        const Real alpha = TBase::residsq_ / normVectorSq(tmp);
        addRescaled(out, dv, alpha);
        addRescaled(rv, tmp2, -alpha);
        const Real residsqnew = normVectorSq(rv);
        const Real beta = residsqnew / TBase::residsq_;
#ifdef DEBUG_INVERTER
        std::cout<<"Testing Operator DAG "<<normVectorSq(tmp)<<" "<<scalarProd(dv,tmp2)<<" "<<residsqnew<<" "<<TBase::residsq_<<" "<<alpha<<" "<<beta<<std::endl;
#endif
        rescaleAdd(dv, beta, rv);
        TBase::residsq_ = residsqnew;
        if (TBase::project_ != 0 && ite % TBase::project_ == 0) {
          projectTracelessVector(out);
          projectTracelessVector(dv);
          projectTracelessVector(rv);
        }
      }
      std::cout << "CGSQ did not converge after "<< TBase::maxit_ <<"steps final resid "<< sqrt(TBase::residsq_) << " " << TBase::residgoal_ << std::endl;
      TBase::iterations_ = TBase::maxit_;
      TBase::fail_ = true;
    }
};

struct RhmcWorkSpace{
    using Vector=BfssDiracVector;
    std::vector<Real> shifts={};
    std::vector<Real> prefact={};
    std::vector<Vector> solutions={};
    Vector phi;
};



template<class OP, class OPDAG>
class CgSqMSolver: public CgnrSolver<OP, OPDAG> {
  public:
    using TBase=CgnrSolver<OP,OPDAG>;
    using Vector=typename TBase::Vector;
    using VectorSpace=std::vector<Vector>;
    using TShiftVector=std::vector<Real>;
  private:
    TShiftVector shift_={};
  public:
    CgSqMSolver(const OP &op, const OPDAG &opdag,const TShiftVector& vect) :
        TBase(op, opdag),shift_(vect) {
    }

    CgSqMSolver(const OP &op, const OPDAG &opdag) :
        TBase(op, opdag){
    }

    void setShiftVector(const TShiftVector& v){
      shift_=v;
    }

    void multiply(const Vector& in, VectorSpace& out){
      const Real normInSq(normVectorSq(in));
      TBase::iterations_ = 0;
      TBase::fail_ = false;


      if (shift_.size() == 0)
        return;
      if (shift_.size() > out.size())
        return;

      TBase::residsq_ =normInSq;

      Vector r(in);
      VectorSpace p(out.size());
      for (size_t i(out.size()); i--;) {
        out[i].resize(in.size());
        out[i].setToZero();
        p[i] = r;
      }

      Vector mp(r);
      Vector tmp(r.size());
      TShiftVector zq(shift_.size() - 1, 1.0);
      TShiftVector ap(shift_.size() - 1, 1.0);

      Real alm1(1.0);
      Real beta(0.0);

      for (unsigned int ite(0); ite < TBase::maxit_; ite++) {
#ifdef DEBUG_INVERTER
        std::cout << "CgSqMSolver:" << ite << " " << sqrt(TBase::residsq_)<< " " << TBase::residgoal_ <<" "<<normInSq<< std::endl;
#endif
        if (TBase::residsq_ != TBase::residsq_) {
                std::cout << "CGSQERR:" << TBase::residsq_ << std::endl;
                TBase::fail_=true;
                TBase::iterations_ = ite;
                return;
        }
        if (TBase::residsq_ < TBase::residgoal_ * TBase::residgoal_ * normInSq ) {
#ifdef DEBUG_INVERTER
        std::cout << "CgSqMSolver:Converged " << ite << " " << sqrt(TBase::residsq_)<< " " << TBase::residgoal_ <<" "<<normInSq<< std::endl;
#endif
          TBase::iterations_ = ite;
          return;
        }
        TBase::op_.multiply(p[0],tmp);
        TBase::opdag_.multiply(tmp,mp);
        Real dpr(normVectorSq(tmp));
        if(shift_[0]!=0.0){
          dpr+=shift_[0]*normVectorSq(p[0]);
        }
        const Real alpha = TBase::residsq_ / dpr;
        addRescaled(out[0], p[0], alpha);
        addRescaled(r, mp, -alpha);
        if(shift_[0]!=0.0){
          addRescaled(r, p[0], -alpha*shift_[0]);
        }
        const Real residsqnew = normVectorSq(r);

        for (size_t s = 0; s < shift_.size() - 1; s++) {
          const Real sh(shift_[s + 1] - shift_[0]);
          const Real invz(
              (1.0 + alpha * sh) + (alpha / alm1) * beta * (1 - zq[s]));
          zq[s] = 1.0 / invz;
          ap[s] = (alpha / alm1) * (zq[s] * ap[s]);
          addRescaled(out[s + 1], p[s + 1], ap[s]);
        }
        beta = residsqnew / TBase::residsq_;
        rescaleAdd(p[0],beta,r);
        for (size_t s = 0; s < shift_.size() - 1; s++) {
          rescaleAdd(p[s + 1], beta * zq[s], r);
        }
        alm1 = alpha;
        TBase::residsq_ = residsqnew;
      }
      std::cout << "CgSqMSolver did not converge " << sqrt(TBase::residsq_) << " " << TBase::residgoal_ << std::endl;
      TBase::iterations_ = TBase::maxit_;
      TBase::fail_ = true;
    }

};

#endif /* SRC_BFSS_CGMSOLVER_H_ */
