/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_BFSSBOSONICACTION_H_
#define SRC_BFSS_BFSSBOSONICACTION_H_

#include "bfsssystemstate.h"
#include "bfssconfig.h"


class BfssBosonicAction {
  private:
    Real kinaction_ = 0.0;
    Real comaction_ = 0.0;
    Real gaugefix_ = 0.0;
    Real alphaconstraint_ = 0.0;
    Real trx2_ = 0.0;
    Real trx2Constraints_ = 0.0;
    Real bmnpart_=0.0;
  public:
    using Mat=BfssConfiguration::MatrixType;
    using State=BfssSystemState;
    using Config=BfssConfiguration;
    BfssBosonicAction() {
    }
    BfssBosonicAction(const BfssBosonicAction &in) = default;
    BfssBosonicAction& operator=(const BfssBosonicAction &in) = default;
    ~BfssBosonicAction() {
    }

    void update(const State &state) {
      const Config &conf(state.config());
      const Real latticeSpacing(state.parameters().latticeSpacing());
      const Real invLatticeSpacing(1.0 / state.parameters().latticeSpacing());
      Real msquare(2.0 * state.parameters().bosonMassSquare() * latticeSpacing * latticeSpacing);
      Real coeff1(-2.0);
      const Real coeff2(3.0/2.0);
      if(state.parameters().improvedaction){
        msquare+=13.0/2.0;
        coeff1=-8.0;
      }else{
        msquare+=2.0;
      }
      Real kinpart(0);
      Real commpart(0);
      Real traceX2(0);
      Real bmnpart1(0);
      Real bmnpart2(0);
#pragma omp parallel for reduction(+:kinpart,commpart,traceX2,bmnpart1,bmnpart2)
      for (size_t site = 0; site < conf.size(); site++) {
        const size_t up = conf.up(site);
        const size_t up2 = conf.up(up);
        Mat tmp,com;
        for (size_t e(Config::elementnum); e--;) {

          const Real trx2(real((conf(site, e) * conf(site, e)).trace()));
          kinpart += msquare * trx2;
          traceX2 += trx2;

            if (state.parameters().ungauged) {
              kinpart += coeff1 * real((conf(up, e) * conf(site, e)).trace());
            } else {
              kinpart += coeff1 * real((conf.expPhases() * conf(up, e) * conf.expConjPhases() * conf(site, e)).trace());
            }

          if(state.parameters().improvedaction){ // improved action
           if (state.parameters().ungauged) {
             kinpart += coeff2 * real((conf(up2, e) * conf(site, e)).trace());
           } else {
             kinpart += coeff2 * real((conf.expPhases(2) * conf(up2, e) * conf.expConjPhases(2) * conf(site, e)).trace());
           }
          } // end improved action

          for (size_t l(e); l--;) {
            //Mat com(conf(site, e) * conf(site, l)-conf(site, l) * conf(site, e));
            // Assume hermiticity.
            tmp.noalias()=conf(site, e) * conf(site, l);
            com.noalias()=tmp-tmp.adjoint();
            commpart -= real((com * com).trace());
            if(state.parameters().mu!=0.0&&e==2&&l==1){ // part one of BMN model
              bmnpart2+=-real(Complex(0.0,3.0)*(conf(site,0)*com).trace());
            }
          }
          if(state.parameters().mu!=0.0){// part two of BMN model
            if(e<3){
              bmnpart1+=trx2/2.0;
            }else{
              bmnpart1+=trx2/8.0;
            }
          }
        }
      }

      trx2_ = traceX2 / Real(numcolor * conf.size());
      kinaction_ = 0.5 * static_cast<Real>(numcolor) * invLatticeSpacing * kinpart*state.parameters().kinPrefact;
      comaction_ = 0.5 * static_cast<Real>(numcolor) * latticeSpacing * commpart*state.parameters().commPrefact;

      if(state.parameters().mu!=0.0){ // BMN part
         bmnpart_=static_cast<Real>(numcolor) * latticeSpacing *(state.parameters().mu*(state.parameters().mu*bmnpart1+bmnpart2));
      }

      // now gauge part
      gaugefix_ = 0.0;
      alphaconstraint_ = 0.0;
      if (!state.parameters().ungauged) {
        for (size_t n(0); n < numcolor; n++)
          for (size_t m(n + 1); m < numcolor; m++) {
            gaugefix_ -= 2.0 * std::log(std::fabs(std::sin(0.5 * (conf.phases()(n) - conf.phases()(m))+state.parameters().gaugefixRegulator)));
          }
        gaugefix_*=state.parameters().alphaFixPrefact;
        const Real maxdiffalpha(conf.phases().maxCoeff() - conf.phases().minCoeff());
        if (maxdiffalpha < 2.0 * pi) {
          alphaconstraint_ = -std::log(2.0 * pi - maxdiffalpha+state.parameters().alphaConstraintRegulator);
        }else{
          // This is not present in the fortran code, but makes it more consistent.
          alphaconstraint_=state.parameters().alphaConstraintRepulsion*(maxdiffalpha-2.0*pi);
        }
        alphaconstraint_*=state.parameters().alphaConstraintPrefact;
      }

      trx2Constraints_ = 0.0;
      if (state.parameters().trX2CutG > 0 && trx2_ > state.parameters().trX2Cut && state.parameters().kinPrefact!=0.0) {
        trx2Constraints_ = state.parameters().trX2CutG * (trx2_ - state.parameters().trX2Cut) * Real(numcolor);
      }
    }

    void addForce(const State &state, Config &force, const Real prefact) {
      const Config &conf(state.config());
      const Real latticeSpacing(state.parameters().latticeSpacing());
      const Real invLatticeSpacing(1.0 / state.parameters().latticeSpacing());
      Real msquare(2.0 * state.parameters().bosonMassSquare() * latticeSpacing * latticeSpacing);
      const Real forcefactKin(prefact*state.parameters().kinPrefact * 0.5 * static_cast<Real>(numcolor) / latticeSpacing);
      const Real forcefactComm(prefact*state.parameters().commPrefact * 0.5 * static_cast<Real>(numcolor) * latticeSpacing);
      const Real bmnprefact1(prefact*static_cast<Real>(numcolor) * latticeSpacing *state.parameters().mu*state.parameters().mu);
      const Real bmnprefact2(prefact*static_cast<Real>(numcolor) * latticeSpacing *state.parameters().mu);

      Real coeff1(-2.0);
      const Real coeff2(3.0/2.0);

      if(state.parameters().improvedaction){
        msquare+=13.0/2.0;
        coeff1=-8.0;
      }else{
        msquare+=2.0;
      }

      if (state.parameters().trX2CutG > 0) {
        const Real trx2(sumTraceX2Norm(conf));
        if (trx2 > state.parameters().trX2Cut&& state.parameters().kinPrefact!=0.0) {
          msquare += state.parameters().trX2CutG / Real(conf.size()) / forcefactKin;
        }
      }

#pragma omp parallel for
      for (size_t site = 0; site < force.size(); site++) {
        const size_t up = force.up(site);
        const size_t dn = force.dn(site);
        const size_t up2 = force.up(up);
        const size_t dn2 = force.dn(dn);
        Mat tmp,com;
        for (size_t e(Config::elementnum); e--;) {

          if(state.parameters().ungauged){
          force(site, e).noalias() += forcefactKin
              * (msquare * conf(site, e)
               +0.5*coeff1*(conf(up, e) +  conf(dn, e)));
          }else{
            force(site, e).noalias() += forcefactKin
                * (msquare * conf(site, e)
                +0.5*coeff1*( conf.expPhases() * conf(up, e) * conf.expConjPhases() + conf.expConjPhases() * conf(dn, e) * conf.expPhases()));
          }

          if(state.parameters().improvedaction){ // improved action
            if(state.parameters().ungauged){
            force(site, e).noalias() += forcefactKin
                * 0.5*coeff2*(conf(up2, e) +  conf(dn2, e));
            }else{
              force(site, e).noalias() += forcefactKin
                  * 0.5*coeff2*( conf.expPhases(2) * conf(up2, e) * conf.expConjPhases(2) + conf.expConjPhases(2) * conf(dn2, e) * conf.expPhases(2));
            }
          } // end improved action

          for (size_t l(Config::elementnum); l--;) {
            // assumes hermiticiy
            tmp.noalias()=conf(site, e) * conf(site, l);
            com.noalias()=conf(site, l)*(tmp-tmp.adjoint());
            force(site, e).noalias() -= forcefactComm * (com + com.adjoint());
            // otherwise:
            //Mat com(conf(site, e) * conf(site, l)-conf(site, l) * conf(site, e));
            //force(site, e).noalias() -= forcefactComm * (conf(site, l) * com - com * conf(site, l));
            if(state.parameters().mu!=0.0){ // BMN part 1
              // hard coded
              if(e==1&&l==2){
                force(site, 0).noalias() +=0.5*bmnprefact2*Complex(0.0,3.0)*(tmp-tmp.adjoint());
              }
              if(e==2&&l==0){
                force(site, 1).noalias() +=0.5*bmnprefact2*Complex(0.0,3.0)*(tmp-tmp.adjoint());
              }
              if(e==0&&l==1){
                force(site, 2).noalias() +=0.5*bmnprefact2*Complex(0.0,3.0)*(tmp-tmp.adjoint());
              }
            }
          }

          if(state.parameters().mu!=0.0){ // BMN part 2
            if(e<3){
              force(site, e).noalias() +=bmnprefact1/2.0*conf(site, e);
            }else{
              force(site, e).noalias() +=bmnprefact1/8.0*conf(site, e);
            }
          }
        }
      }

      // gauge force part
      if (!state.parameters().ungauged && !state.parameters().backroundgauge) {
        GaugeVector tmp;
        tmp.setZero();
        for (size_t site = 0; site < force.size(); site++) {
          const size_t up = force.up(site);
          const size_t up2 = force.up(up);
          for (size_t e(Config::elementnum); e--;) {
            tmp+=0.5*coeff1*(conf.expPhases()*conf(up, e)  * conf.expConjPhases() *conf(site, e)).diagonal();
            if(state.parameters().improvedaction){
              tmp+=0.5*coeff2*2.0*(conf.expPhases(2)*conf(up2, e)  * conf.expConjPhases(2) *conf(site, e)).diagonal();
            }
          }
        }
        const Complex isize(0.0,1.0/static_cast<Real>(force.size()));
        for (size_t n(numcolor); n--;) {
          force.phases()(n) += forcefactKin * 2.0*real(isize *tmp(n));
        }
        for (size_t n(0); n < numcolor; n++)
          for (size_t m(0); m < numcolor; m++)
            if (m != n) {
              force.phases()(n) += -0.5 / (std::tan(0.5 * (conf.phases()(n) - conf.phases()(m))+(m>n ? 1.0:-1.0)*state.parameters().gaugefixRegulator))*state.parameters().alphaFixPrefact;
            }

        auto maxInd=maxIndex(conf.phases());
        auto minInd=minIndex(conf.phases());
        const Real maxdiffalpha(maxInd.first-minInd.first);
        if (maxdiffalpha < 2.0 * pi) {
          force.phases()(maxInd.second) += (0.5/(2.0 * pi - maxdiffalpha+state.parameters().alphaConstraintRegulator))*state.parameters().alphaConstraintPrefact;
          force.phases()(minInd.second) += (-0.5/(2.0 * pi - maxdiffalpha+state.parameters().alphaConstraintRegulator))*state.parameters().alphaConstraintPrefact;
        }else{
          force.phases()(maxInd.second)+=0.5*state.parameters().alphaConstraintRepulsion;
          force.phases()(minInd.second)-=0.5*state.parameters().alphaConstraintRepulsion;
        }
      }

    }

    Real action() const {
      return (kinaction_ + comaction_ + gaugefix_ + alphaconstraint_ + trx2Constraints_+bmnpart_);
    }
    Real kineticPart() const {
      return (kinaction_);
    }
    Real commutatorPart() const {
      return (comaction_);
    }

    Real gaugeFix() const {
      return gaugefix_;
    }

    Real bmnPart() const {
      return bmnpart_;
    }

    void logActionParts(std::ostream& os)const{
      os<<"BosonicAction:"<<action()<<":kin:"<<kineticPart()<<":comm:"<<commutatorPart()<<":bmn:"<<bmnPart()<<":gf:"<<gaugeFix()<<":Aconstr:"<<alphaconstraint_<<":x2const:"<<trx2Constraints_<<":x2:"<<trx2_;
    }

};

#endif /* SRC_BFSS_BFSSBOSONICACTION_H_ */
