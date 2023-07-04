/*
 * diracoperator.h
 *
 *  Created on: 27.03.2023
 *      Author: Georg Bergner
 */

#ifndef DIRACOPERATOR_H_
#define DIRACOPERATOR_H_

#include "bfssconfig.h"

struct BfssDiracOperatorParameter {
    // additional input:
    Real mu=0.0;
    Real fermmass=0.0;
    bool improved=false;
    bool gauged=true;
    // this is not the most efficient implementation of fermion boundary conditions, but it is not important for optimizations at the moment.
    Real fermionbc=1.0;

    // set by functions below.
    Complex prefactMu=1.0;
    Complex prefactorKin = 1.0;
    Complex prefactComm = 1.0;
    void setDefaults(const Real latticespacing) {
      setMasanorisConventions(latticespacing);
    }
    void setTextBook(const Real latticespacing) {
      prefactorKin = Complex(0.0, Real(numcolor));
      prefactComm = -latticespacing * Complex(Real(numcolor), 0.0);
      prefactMu=Complex(0.0,0.75*mu*Real(numcolor)*latticespacing);
    }
    void setMasanorisConventions(const Real latticespacing) {
      prefactorKin = Complex(1.0, 0.0);
      prefactComm = -latticespacing * Complex(0.0, -1.0);
      prefactMu=Complex(-0.75*mu*latticespacing,0.0);
    }
};

template<unsigned int mu>
struct SigmaPauli;

template<>
struct SigmaPauli<0> {
    static size_t index(size_t ind) {
      return ind == 0 ? 1 : 0;
    }
    static Complex value(const size_t ind) {
      return 1.0;
    }
};

template<>
struct SigmaPauli<1> {
    static size_t index(size_t ind) {
      return ind == 0 ? 1 : 0;
    }
    static Complex value(const size_t ind) {
      return Complex(0.0, (ind == 0 ? -1.0 : 1.0));
    }
};

template<>
struct SigmaPauli<2> {
    static size_t index(size_t ind) {
      return ind;
    }
    static Complex value(const size_t ind) {
      return (ind == 0 ? 1.0 : -1.0);
    }
};

// Sigma 4 is defined as identity
template<>
struct SigmaPauli<4> {
    static size_t index(size_t ind) {
      return ind;
    }
    static Complex value(const size_t ind) {
      return 1.0;
    }
};

// Sigma 5 is defined as -i*identity
template<>
struct SigmaPauli<5> {
    static size_t index(size_t ind) {
      return ind;
    }
    static Complex value(const size_t ind) {
      return Complex(0.0, -1.0);
    }
};

template<unsigned int s1, unsigned int s2, unsigned int s3, unsigned int s4>
struct SigmaTensor {
    static size_t index(const size_t ind) {
      const size_t ind1 = ind / 8;
      const size_t ind2 = (ind - ind1 * 8) / 4;
      const size_t ind3 = (ind - ind1 * 8 - ind2 * 4) / 2;
      const size_t ind4 = (ind - ind1 * 8 - ind2 * 4 - ind3 * 2);
      return SigmaPauli<s1>::index(ind1) * 8 + SigmaPauli<s2>::index(ind2) * 4 + SigmaPauli<s3>::index(ind3) * 2 + SigmaPauli<s4>::index(ind4);
    }
    static Complex value(const size_t ind) {
      const size_t ind1 = ind / 8;
      const size_t ind2 = (ind - ind1 * 8) / 4;
      const size_t ind3 = (ind - ind1 * 8 - ind2 * 4) / 2;
      const size_t ind4 = (ind - ind1 * 8 - ind2 * 4 - ind3 * 2);
      return SigmaPauli<s1>::value(ind1) * SigmaPauli<s2>::value(ind2) * SigmaPauli<s3>::value(ind3) * SigmaPauli<s4>::value(ind4);
    }
    static Complex valueAdjoint(const size_t ind) {
      return value(ind);
    }
};

template<unsigned int mu>
struct Gamma10D;

// gamma1=sigma_3 x id x id x id
template<>
struct Gamma10D<0> : SigmaTensor<2, 4, 4, 4> {
};

// gamma2=sigma_2 x sigma_2 x sigma_2 x sigma_2
template<>
struct Gamma10D<1> : SigmaTensor<1, 1, 1, 1> {
};

template<>
struct Gamma10D<2> : SigmaTensor<1, 1, 4, 0> {
};

template<>
struct Gamma10D<3> : SigmaTensor<1, 1, 4, 2> {
};

template<>
struct Gamma10D<4> : SigmaTensor<1, 0, 1, 4> {
};

template<>
struct Gamma10D<5> : SigmaTensor<1, 2, 1, 4> {
};

template<>
struct Gamma10D<6> : SigmaTensor<1, 4, 0, 1> {
};

template<>
struct Gamma10D<7> : SigmaTensor<1, 4, 2, 1> {
};

template<>
struct Gamma10D<8> : SigmaTensor<5, 4, 4, 4> {
    static Complex valueAdjoint(const size_t ind) {
       return -SigmaTensor<5, 4, 4, 4>::value(ind);
     }
};

template<>
struct Gamma10D<9> : SigmaTensor<0, 4, 4, 4> {

};

struct Gamma123 {
    static size_t index(const size_t ind) {
      return Gamma10D<0>::index(Gamma10D<1>::index(Gamma10D<2>::index(ind)));
    }
    static Complex value(const size_t ind) {
      return Gamma10D<0>::value(Gamma10D<1>::index(Gamma10D<2>::index(ind))) * Gamma10D<1>::value(Gamma10D<2>::index(ind)) * Gamma10D<2>::value(ind);
    }
};

class BfssDiracOperator {
  public:
    using Config=BfssConfiguration;
    using Vector=BfssDiracVector;
  private:
    using Mat=BfssConfiguration::MatrixType;
    bool dag=false;

    const BfssDiracOperatorParameter *parameter_ = nullptr;
    const Config *config_ = nullptr;
    const Config& conf() {
      return *config_;
    }
    const BfssDiracOperatorParameter parameters() {
      return *parameter_;
    }
  public:
    BfssDiracOperator(const Config *config, const BfssDiracOperatorParameter *parameter) :
        parameter_(parameter), config_(config) {

    }

    BfssDiracOperator(const BfssDiracOperator& in)=default;

    BfssDiracOperator& operator=(const BfssDiracOperator& in)=default;

    void setDag(const bool v){
      dag=v;
    }

    // be careful= always assuming no alias (don't use multiply(in,in)).
    void multiply(const Vector &in, Vector &out) {
      multiplyAllKinetic(in,out);
      multiplyCommutatorPart(in, out);
      if(parameters().mu!=0.0){
        multiplyBmnPart(in,out);
      }
      if(parameters().fermmass!=0.0){
        multiplyMassTerm(in,out);
      }
    }

    void addForce(Config& in,const Vector& X, const Vector& Y, const Real pref){
      if(parameters().gauged){
       if(parameters().improved){
         addForceKineticImproved(in,X,Y,pref);
       }else{
         addForceKinetic(in,X,Y,pref);
       }
      }
      addForceCommutatorPart(in,X,Y,pref);
      // no force contribution from BMN and mass part.
    }

  private:

    // not the most efficient implementation, but not crucial at the moment.
    Real boundaryFactorUp(const size_t site){
      return site+1==conf().size() ? parameters().fermionbc:1.0;
    }

    Real boundaryFactorDown(const size_t site){
      return site==0 ? parameters().fermionbc:1.0;
    }

    void multiplyAllKinetic(const Vector &in, Vector &out){
      if(parameters().gauged){
          multiplyKineticImproved(in,out);
        if(parameters().improved){
        }else{
          multiplyKinetic(in, out);
        }
      }else{ // ungauged
        if(parameters().improved){
          multiplyKineticImprovedUngauged(in,out);
        }else{
          multiplyKineticUngauged(in, out);
        }
      }
    }

    void multiplyKineticUngauged(const Vector &in, Vector &out) {
      const Complex prefact(dag ? -conj(parameters().prefactorKin):parameters().prefactorKin);
#pragma omp parallel for
      for (size_t site = 0; site < in.size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        // Hard coded gamma_10
        // Please be aware that D_+ and D_- are exchanged compared to the notes.
        // This is the same as in the fortran code. (same for improved version).
        for (size_t al(0); al < 8; al++) {
          out(site, al + 8).noalias() = prefact * (bcfactup*in(up, al)  - in(site, al)); //D_+
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        for (size_t al(0); al < 8; al++) {
          out(site, al).noalias() = prefact * (in(site, al + 8) - bcfactdn*in(dn, al + 8)); //D_-
        }
      }
    }

    void multiplyKinetic(const Vector &in, Vector &out) {
      const Complex prefact(dag ? -conj(parameters().prefactorKin):parameters().prefactorKin);
#pragma omp parallel for
      for (size_t site = 0; site < in.size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        // Hard coded gamma_10
        // Please be aware that D_+ and D_- are exchanged compared to the notes.
        // This is the same as in the fortran code. (same for improved version).
        for (size_t al(0); al < 8; al++) {
          out(site, al + 8).noalias() = prefact * (bcfactup*conf().expPhases() * in(up, al) * conf().expConjPhases() - in(site, al)); //D_+
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        for (size_t al(0); al < 8; al++) {
          out(site, al).noalias() = prefact * (in(site, al + 8) - bcfactdn*conf().expConjPhases() * in(dn, al + 8) * conf().expPhases()); //D_-
        }
      }
    }

    void addForceKinetic(Config& in,const Vector& X, const Vector& Y, const Real pref){
      const Complex prefact(pref*parameters().prefactorKin);
      GaugeVector tmp;
      tmp.setZero();
      for (size_t site = 0; site < conf().size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        for (size_t al(0); al < 8; al++) {
          tmp+=(bcfactup*prefact)*(conf().expPhases()*X(up, al)  * conf().expConjPhases()*Y(site, al + 8).adjoint()).diagonal();
          tmp-=(bcfactup*prefact)*(conf().expConjPhases()*Y(site, al + 8).adjoint()*conf().expPhases()*X(up, al)).diagonal();
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        for (size_t al(0); al < 8; al++) {
          tmp+=(bcfactdn*prefact)*(conf().expConjPhases()*X(dn, al + 8)  * conf().expPhases()*Y(site, al).adjoint()).diagonal();
          tmp-=(bcfactdn*prefact)*(conf().expPhases()*Y(site, al).adjoint()*conf().expConjPhases()*X(dn, al + 8)).diagonal();
        }
      }
      const Complex isize(0.0,1.0/static_cast<Real>(in.size()));
      for (size_t n(numcolor); n--;) {
        in.phases()(n) += -real(isize *tmp(n)); // takes into account D^dag contribution.
      }

    }

    void multiplyKineticImprovedUngauged(const Vector &in, Vector &out) {
      const Complex prefact(dag ? -conj(parameters().prefactorKin):parameters().prefactorKin);
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        const size_t up2 = conf().up(up);
        const Real bcfactup2=boundaryFactorUp(up)*bcfactup;
        // Hard coded gamma_10
        for (size_t al(0); al < 8; al++) {
          out(site, al + 8).noalias() = prefact * (
              -(0.5*bcfactup2)*in(up2, al)
              +(2.0*bcfactup)*in(up, al)
              -1.5*in(site, al)); //D_+
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        const size_t dn2 = conf().dn(dn);
        const Real bcfactdn2=boundaryFactorDown(dn)*bcfactdn;
        for (size_t al(0); al < 8; al++) {
          out(site, al).noalias() = prefact * (
              1.5*in(site, al + 8)
              -(2.0*bcfactdn)*in(dn, al + 8)
              +(0.5 *bcfactdn2)* in(dn2, al + 8)); //D_-
        }
      }
    }

    void multiplyKineticImproved(const Vector &in, Vector &out) {
      const Complex prefact(dag ? -conj(parameters().prefactorKin):parameters().prefactorKin);
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        const size_t up2 = conf().up(up);
        const Real bcfactup2=boundaryFactorUp(up)*bcfactup;
        // Hard coded gamma_10
        for (size_t al(0); al < 8; al++) {
          out(site, al + 8).noalias() = prefact * (
              -(0.5*bcfactup2)*conf().expPhases(2) * in(up2, al) * conf().expConjPhases(2)
              +(2.0*bcfactup)*conf().expPhases(1) * in(up, al) * conf().expConjPhases(1)
              -1.5*in(site, al)); //D_+
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        const size_t dn2 = conf().dn(dn);
        const Real bcfactdn2=boundaryFactorDown(dn)*bcfactdn;
        for (size_t al(0); al < 8; al++) {
          out(site, al).noalias() = prefact * (
              1.5*in(site, al + 8)
              -(2.0*bcfactdn)*conf().expConjPhases(1) * in(dn, al + 8) * conf().expPhases(1)
              +(0.5 *bcfactdn2)*conf().expConjPhases(2) * in(dn2, al + 8) * conf().expPhases(2)); //D_-
        }
      }
    }

    void addForceKineticImproved(Config& in,const Vector& X, const Vector& Y, const Real pref){
      const Complex prefact(pref*parameters().prefactorKin);
      GaugeVector tmp;
      tmp.setZero();
      for (size_t site = 0; site < conf().size(); site++) {
        const size_t up = conf().up(site);
        const Real bcfactup=boundaryFactorUp(site);
        const size_t up2 = conf().up(up);
        const Real bcfactup2=boundaryFactorUp(up)*bcfactup;
        for (size_t al(0); al < 8; al++) {
          tmp+=-prefact*(0.5*bcfactup2)*2.0*(conf().expPhases(2)*X(up2, al)  * conf().expConjPhases(2)*Y(site, al + 8).adjoint()).diagonal();
          tmp-=-prefact*(0.5*bcfactup2)*2.0*(conf().expConjPhases(2)*Y(site, al + 8).adjoint()*conf().expPhases(2)*X(up2, al)).diagonal();
          tmp+=prefact*(2.0*bcfactup)*(conf().expPhases()*X(up, al)  * conf().expConjPhases()*Y(site, al + 8).adjoint()).diagonal();
          tmp-=prefact*(2.0*bcfactup)*(conf().expConjPhases()*Y(site, al + 8).adjoint()*conf().expPhases()*X(up, al)).diagonal();
        }
        const size_t dn = conf().dn(site);
        const Real bcfactdn=boundaryFactorDown(site);
        const size_t dn2 = conf().dn(dn);
        const Real bcfactdn2=boundaryFactorDown(dn)*bcfactdn;
        for (size_t al(0); al < 8; al++) {
          tmp+=prefact*(2.0*bcfactdn)*(conf().expConjPhases()*X(dn, al + 8)  * conf().expPhases()*Y(site, al).adjoint()).diagonal();
          tmp-=prefact*(2.0*bcfactdn)*(conf().expPhases()*Y(site, al).adjoint()*conf().expConjPhases()*X(dn, al + 8)).diagonal();
          tmp+=-prefact*(0.5 *bcfactdn2)*2.0*(conf().expConjPhases(2)*X(dn2, al + 8)  * conf().expPhases(2)*Y(site, al).adjoint()).diagonal();
          tmp-=-prefact*(0.5 *bcfactdn2)*2.0*(conf().expPhases(2)*Y(site, al).adjoint()*conf().expConjPhases(2)*X(dn2, al + 8)).diagonal();
        }
      }
      const Complex isize(0.0,1.0/static_cast<Real>(in.size()));
      for (size_t n(numcolor); n--;) {
        in.phases()(n) += -real(isize *tmp(n));
      }

    }

    template<int n, class V>
    struct LoopMu {
        static void exec(const size_t site, const V &in, V &out, const Complex &fact, const Config &config,const bool dag) {
          for (size_t al(0); al < V::elementnum; al++) {
            const Complex gval(dag ? Gamma10D<n>::valueAdjoint(al):Gamma10D<n>::value(al));
            const Complex val((dag ? conj(fact):fact) * gval);
            const size_t ind(Gamma10D<n>::index(al));
            out(site, al).noalias() += val * (config(site, n) * in(site, ind) - in(site, ind) * config(site, n));
          }
          LoopMu<n - 1, V>::exec(site, in, out, fact, config,dag);
        }
        static void addForce(const size_t site,Config& out,const Vector& X,const Vector& Y, const Complex& fact,const Config& config,const bool dag){
          Mat tmp;
          for (size_t al(0); al < V::elementnum; al++) {
            const Complex gval(dag ? Gamma10D<n>::valueAdjoint(al):Gamma10D<n>::value(al));
            const Complex val((dag ? conj(fact):fact) * gval);
            const size_t ind(Gamma10D<n>::index(al));
            tmp.noalias()=val *(X(site,ind)*Y(site, al).adjoint()-Y(site, al).adjoint()*X(site,ind));
            out(site, n).noalias() +=-0.5*(tmp +tmp.adjoint()); // takes into account D^dag contribution.
          }
          LoopMu<n - 1, V>::addForce(site,out,X,Y, fact, config,dag);
        }
    };
    template<class V>
    struct LoopMu<-1, V> {
        static void exec(const size_t site, const V &in, V &out, const Complex &fact, const Config &config,const bool dag) {
        }
        static void addForce(const size_t site,Config& out,const Vector& X,const Vector& Y, const Complex& fact,const Config& config,const bool dag){
        }
    };

    void multiplyMassTerm(const Vector &in, Vector &out) {
      const Complex prefact(parameters().prefactComm);
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        for (size_t al(0); al < Vector::elementnum; al++) {
         out(site, al).noalias() += parameters().fermmass*in(site, al);//conf()(site, 0) * in(site, al)
        }
      }
    }

/*    void addFakePartForce(Config& out,const Vector& X, const Vector& Y, const Real pref){
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        for (size_t al(0); al < Vector::elementnum; al++) {
//         out(site, 0).noalias() += -0.5*pref*(X(site, al)*Y(site,al).adjoint()+Y(site, al)*X(site,al).adjoint());
        }
      }
    }*/

    void multiplyCommutatorPart(const Vector &in, Vector &out) {
      const Complex prefact(parameters().prefactComm);
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        LoopMu<Config::elementnum - 1, Vector>::exec(site, in, out, prefact, conf(),dag);
      }
    }

    void addForceCommutatorPart(Config& out,const Vector& X, const Vector& Y, const Real pref){
      const Complex prefact(pref*parameters().prefactComm);
#pragma omp parallel for
      for (size_t site = 0; site < conf().size(); site++) {
        LoopMu<Config::elementnum - 1, Vector>::addForce(site, out, X,Y, prefact, conf(),dag);
      }
    }

    void multiplyBmnPart(const Vector &in, Vector &out) {
      const Complex prefact(dag ? -conj(parameters().prefactMu):parameters().prefactMu);
#pragma omp parallel for
      for (size_t site = 0; site < in.size(); site++) {
        for (size_t al(0); al < Vector::elementnum; al++) {
          out(site,al).noalias()+=prefact*Gamma123::value(al)*in(site,Gamma123::index(al));
        }
      }
    }

};

#endif /* DIRACOPERATOR_H_ */
