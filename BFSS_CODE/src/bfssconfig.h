/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef BFSSCONFIG_H_
#define BFSSCONFIG_H_
#include "basicdef.h"
#include "RandomNumberGenerator.h"
#include <fstream>
#include <list>
#include <iomanip>
#ifdef DYNAMIC_ALLOCATION
#define USE_DYNAMIC_ARRAY
#endif

template<unsigned int ESIZE>
class MatrixModelBaseConfiguration {
public:
  static const size_t elementnum = ESIZE;
  static const size_t elementsize = numcolor * numcolor * elementnum;
#ifndef USE_DYNAMIC_ARRAY
  typedef GaugeMatrix SiteElement[elementnum];
  typedef GaugeMatrix StoredGaugeMatrix;
  typedef GaugeMatrix& ReturnType;
  typedef GaugeMatrix MatrixType;
  typedef const GaugeMatrix& ConstReturnType;
  ReturnType mapElement(StoredGaugeMatrix& in) {
    return in;
  }
  ConstReturnType mapElement(const StoredGaugeMatrix& in) const {
    return in;
  }
#else
  typedef DynamicGaugeMatrix MatrixType;
  typedef Complex StoredGaugeMatrix[numcolor*numcolor];
  typedef StoredGaugeMatrix SiteElement[elementnum];
  typedef GaugeMatrixMap ReturnType;
  typedef GaugeMatrixMapConst ConstReturnType;
  ReturnType mapElement(StoredGaugeMatrix& in) {
    return (ReturnType(in,numcolor,numcolor));
  }
  ConstReturnType mapElement(const StoredGaugeMatrix& in)const {
    return (ConstReturnType(in,numcolor,numcolor));
  }
#endif

private:
  SiteElement* data_=nullptr;
  size_t size_=0;
  void allocate_(){
    data_ = new SiteElement[size_];
  }
  void free_(){
    delete[] data_;
  }
  void copyData_(const MatrixModelBaseConfiguration& in){
#pragma omp parallel for
    for (size_t i=0;i<size(); i++) {
      for (size_t m(elementnum); m--;)
        this->operator()(i,m) = in(i,m);
    }
  }
  void resize_(const size_t size){
    if (size_!=size) {
      free_();
      size_=size;
      allocate_();
    }
  }
public:

  size_t size()const{
    return size_;
  }
  size_t up(const size_t site)const{
    return (site+1)%size_;
  }

  size_t dn(const size_t site)const{
    return (site+(size_-1))%size_;
  }

  MatrixModelBaseConfiguration() {
  }

  MatrixModelBaseConfiguration(const size_t size) {
    size_=size;
   allocate_();
  }


  virtual ~MatrixModelBaseConfiguration() {
    free_();
  }

  void resize(const size_t size){
    if(size!=size_){
     size_=size;
     free_();
     allocate_();
    }
  }



  MatrixModelBaseConfiguration(const MatrixModelBaseConfiguration& in):size_(in.size_) {
    allocate_();
    copyData_(in);
  }

  const MatrixModelBaseConfiguration& operator=(
      const MatrixModelBaseConfiguration& in) {
    resize_(in.size());
    copyData_(in);
    return (*this);
  }


  ReturnType operator()(const size_t site, const size_t c) {
    return (mapElement(data_[site][c]));
  }

  ConstReturnType operator()(const size_t site, const size_t c) const {
    return (mapElement(data_[site][c]));
  }

  void setToIdentity() {
#pragma omp parallel for
    for (size_t i=0;i<size(); i++) {
      for (size_t mu(elementnum); mu--;) {
        mapElement(data_[i][mu]).setIdentity();
      }
    }
  }

  void setToZero() {
#pragma omp parallel for
    for (size_t i=0;i<size(); i++) {
      for (size_t mu(elementnum); mu--;) {
        mapElement(data_[i][mu]).setZero();
      }
    }
  }


};

template<unsigned int FILEDNUM>
class MatrixModelConfiguration: public MatrixModelBaseConfiguration<FILEDNUM> {
public:
    using TBase=MatrixModelBaseConfiguration<FILEDNUM>;
private:
  RealGaugeVector phases_={};
  void copyDataPhase_(const MatrixModelConfiguration& in){
    phases_=in.phases_;
  }

public:


  MatrixModelConfiguration():TBase(){
  }

  MatrixModelConfiguration(const size_t size):TBase(size) {
  }


  ~MatrixModelConfiguration() {

  }


  MatrixModelConfiguration(const MatrixModelConfiguration& in):TBase(in) {
    copyDataPhase_(in);
  }

  const MatrixModelConfiguration& operator=(
      const MatrixModelConfiguration& in) {
    TBase::operator=(in);
    copyDataPhase_(in);
    return (*this);
  }


  auto& phases() {
    return phases_;
  }

  const auto& phases() const{
    return phases_;
  }

  auto expPhases(const int n=1)const{
    GaugeVector tmp(Eigen::exp((Complex(0.0,Real(n)/static_cast<Real>(TBase::size()))*phases_).array()));
    return DiagonalGaugeMatrix(tmp.asDiagonal());
  }

  auto expConjPhases(const int n=1)const{
    GaugeVector tmp(Eigen::exp((Complex(0.0,-Real(n)/static_cast<Real>(TBase::size()))*phases_).array()));
    return DiagonalGaugeMatrix(tmp.asDiagonal());
  }


  void setPhaseToIdentity() {
      for (size_t i(numcolor); i--;) {
        phases_(i)=0.0;

      }
  }


  void setToZero() {
    TBase::setToZero();
    setPhaseToIdentity();
  }


};

using BfssConfiguration=MatrixModelConfiguration<SCALARFIELDS>;

using BfssDiracVector=MatrixModelBaseConfiguration<BFSSSPINSIZE>;

inline void setGaussianRandomX(BfssConfiguration& in,RandomNumberGenerator& rnd) {
  for (size_t i = 0; i < in.size(); i++)
    for (size_t e(BfssConfiguration::elementnum); e--;) {
      std::vector<Real> randnum(numcolor*numcolor,0);
      drawGaussianRandom(randnum,rnd);
      in(i, e).setZero();
      size_t ind(0);
      for (size_t n(0); n < numcolor; n++){
        for (size_t m(n+1); m < numcolor; m++) {
          const Real r1(randnum[ind++]/2.0);
          const Real r2(randnum[ind++]/2.0);
          const Complex v(r1,r2);
         in(i,e)(n,m)=v;
         in(i,e)(m,n)=std::conj(v);
      }
        in(i,e)(n,n)=randnum[ind++]/std::sqrt(2.0);
      }
    }
}

inline void setGaussianRandomVector(BfssDiracVector& in,RandomNumberGenerator& rnd) {
  for (size_t i = 0; i < in.size(); i++)
    for (size_t e(BfssDiracVector::elementnum); e--;) {
      std::vector<Real> randnum(2*numcolor*numcolor,0);
      drawGaussianRandom(randnum,rnd);
      in(i, e).setZero();
      size_t ind(0);
      for (size_t n(0); n < numcolor; n++){
        for (size_t m(0); m < numcolor; m++) {
          const Real r1(randnum[ind++]/2.0);
          const Real r2(randnum[ind++]/2.0);
          const Complex v(r1,r2);
          in(i,e)(n,m)=v;
      }
      }
      const Complex tr1(in(i,e).trace()/Real(numcolor));
      for (size_t n(0); n < numcolor; n++){
      in(i,e)(n,n)-=tr1;
      }
    }
}


inline void setU1RandomVector(BfssDiracVector& in,RandomNumberGenerator& rnd) {
  for (size_t i = 0; i < in.size(); i++)
    for (size_t e(BfssDiracVector::elementnum); e--;) {
      std::vector<Real> randnum(2*numcolor*numcolor,0);
      drawU1Random(randnum,rnd);
      in(i, e).setZero();
      size_t ind(0);
      for (size_t n(0); n < numcolor; n++){
        for (size_t m(0); m < numcolor; m++) {
          const Real r1(randnum[ind++]);
          const Real r2(randnum[ind++]);
          const Complex v(r1,r2);
          in(i,e)(n,m)=v;
      }
      }
      const Complex tr1(in(i,e).trace()/Real(numcolor));
      for (size_t n(0); n < numcolor; n++){
      in(i,e)(n,n)-=tr1;
      }
    }
}


inline void projectTracelessVector(BfssDiracVector& in) {
#pragma omp parallel for
  for (size_t i = 0; i < in.size(); i++)
    for (size_t e(BfssDiracVector::elementnum); e--;) {
      const Complex tr1(in(i,e).trace()/Real(numcolor));
      for (size_t n(0); n < numcolor; n++){
        in(i,e)(n,n)-=tr1;
      }
    }
}

inline void sumTraceX2(const BfssConfiguration& conf,std::vector<Real>& traces){
     traces.resize(BfssConfiguration::elementnum);
     for (size_t e(BfssConfiguration::elementnum); e--;) {
       Real traceX(0);
#pragma omp parallel for reduction(+:traceX)
     for (size_t site = 0; site < conf.size(); site++) {
       traceX+=real((conf(site, e)*conf(site, e)).trace());
       }
      traces[e]=traceX;
     }
     }

inline void traceOutAlpha(BfssConfiguration& in) {
   in.phases()=in.phases().array()-(in.phases().sum()/Real(numcolor));
}

inline void reduceToFitConstraint(BfssConfiguration& in) {
  const Real maxdiffalpha(in.phases().maxCoeff() - in.phases().minCoeff());
  if(maxdiffalpha>(2.0*pi-1e-4)){
    in.phases()*=(2.0*pi-1e-4)/maxdiffalpha;
  }
  std::cout<<(in.phases().maxCoeff() - in.phases().minCoeff())/2.0/pi<<std::endl;
}

inline bool constraintViolationAlpha(const BfssConfiguration& in) {
  const Real maxdiffalpha(in.phases().maxCoeff() - in.phases().minCoeff());
  return (maxdiffalpha >= 2.0 * pi) ;
}

// could remove imag part.
inline void traceOutX(BfssConfiguration& conf){
  for (size_t e(BfssConfiguration::elementnum); e--;) {
    Real traceXr(0.0);
    Real traceXi(0.0);
#pragma omp parallel for reduction(+:traceXr,traceXi)
  for (size_t site = 0; site < conf.size(); site++) {
    const Complex tmp=conf(site, e).trace();
    traceXr+=tmp.real();
    traceXi+=tmp.imag();
    }
   const Complex traceX(Complex(traceXr,traceXi)/Real(conf.size()*numcolor));
#pragma omp parallel for
  for (size_t site = 0; site < conf.size(); site++) {
    for(size_t n(0);n<numcolor;n++){
      conf(site,e)(n,n)-=traceX;
    }
  }
  }
}

inline void setGaussianRandomAlpha(BfssConfiguration& in,RandomNumberGenerator& rnd) {
  std::vector<Real> randnum(numcolor,0);
  drawGaussianRandom(randnum,rnd);
  for(size_t n(0);n<numcolor;n++){
    in.phases()(n)=randnum[n]/std::sqrt(2.0);
  }
}

inline void setUniformRandomAlpha(BfssConfiguration& in,RandomNumberGenerator& rnd) {
  for(size_t n(0);n<numcolor;n++){
    in.phases()(n)=rnd()*2.0*pi-pi;
  }
}

inline void setFortranCodeAlpha(BfssConfiguration& in) {
  for(size_t n(0);n<numcolor;n++){
    in.phases()(n)=-1.0+2.0*Real(n+1)/Real(numcolor);
  }
/*  do imat=1,nmat
     alpha(imat)=-1d0+2d0*dble(imat)/dble(nmat)
  end do */
}

// inline void multiplyByConst(BfssConfiguration& inout, Real fact){
//   for (size_t site=0; site< inout.size(); site++){
//     for (size_t e(BfssConfiguration::elementnum); e--;){
//       inout(site,e) = fact*inout(site,e);
//     }
//   }
// }




inline Real sumAlpha2(const BfssConfiguration& in) {
  Real sum(0.0);
  for(size_t n(0);n<numcolor;n++){
    sum+=in.phases()(n)*in.phases()(n);
  }
  return sum;
}

inline Real sumTraceX(const BfssConfiguration& conf){
     Real traceX(0);
#pragma omp parallel for reduction(+:traceX)
     for (size_t site = 0; site < conf.size(); site++) {
       for (size_t e(BfssConfiguration::elementnum); e--;) {
       traceX+=real(conf(site, e).trace());
       }
     }
     return traceX;
}

inline Real sumTraceX2(const BfssConfiguration& conf){
     Real traceX2(0);
#pragma omp parallel for reduction(+:traceX2)
     for (size_t site = 0; site < conf.size(); site++) {
       for (size_t e(BfssConfiguration::elementnum); e--;) {
       traceX2+=real((conf(site, e) * conf(site, e)).trace());
       }
     }
     return traceX2;
}

inline Real normVectorSq(const BfssDiracVector& v){
     Real traceX2(0);
#pragma omp parallel for reduction(+:traceX2)
     for (size_t site = 0; site < v.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
       traceX2+=real((v(site, e) * v(site, e).adjoint()).trace());
       }
     }
     return traceX2;
}

inline Real sumAbsTraceVector(const BfssDiracVector& v){
     Real tracev(0);
#pragma omp parallel for reduction(+:tracev)
     for (size_t site = 0; site < v.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
       tracev+=std::abs(v(site, e).trace());
       }
     }
     return tracev;
}

inline Complex scalarProd(const BfssDiracVector& v1, const BfssDiracVector& v2){
     Real sre(0);
     Real sim(0);
#pragma omp parallel for reduction(+:sre,sim)
     for (size_t site = 0; site < v1.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
       const Complex tmp((v1(site, e).adjoint() * v2(site, e)).trace());
       sre+=tmp.real();
       sim+=tmp.imag();
       }
     }
     return Complex(sre,sim);
}

inline void rescaleAdd(BfssDiracVector& v1,const Real& fact, const BfssDiracVector& v2){
#pragma omp parallel for
     for (size_t site = 0; site < v1.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         v1(site,e)=fact*v1(site,e)+v2(site,e);
       }
     }
}

inline void addRescaled(BfssDiracVector& v1, const BfssDiracVector& v2,const Real& fact){
#pragma omp parallel for
     for (size_t site = 0; site < v1.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         v1(site,e)+=fact*v2(site,e);
       }
     }
}


// normalized according to previous conventions.
inline Real sumTraceX2Norm(const BfssConfiguration& conf){
     return sumTraceX2(conf)/Real(numcolor*conf.size());;
}

inline Real polyakovLoopAlpha(const BfssConfiguration& conf){
 GaugeVector tmp(Eigen::exp((Complex(0.0,1.0)*conf.phases()).array()));
 return std::abs(tmp.sum())/Real(numcolor);
}

inline void addRescaled(BfssConfiguration& inout, const BfssConfiguration& in, const Real fx, const Real falpha){
#pragma omp parallel for
    for (size_t i = 0; i < inout.size(); i++) {
      for (size_t e(BfssConfiguration::elementnum); e--;) {
        inout(i, e) += fx
            * in(i, e);
      }
    }
    for(size_t n(0);n<numcolor;n++){
      inout.phases()(n)+=falpha*in.phases()(n);
    }
}

inline Real sumProductTraceX(const BfssConfiguration& conf1, const BfssConfiguration& conf2){
Real ret(0.0);
#pragma omp parallel for reduction(+:ret)
for (size_t i=0;i<conf1.size(); i++) {
  for (size_t e(BfssConfiguration::elementnum); e--;) {
    ret +=real((conf1(i, e)
                 * conf2(i, e)).trace());
  }
}
return ret;
}

inline Real sumProductAlpha(const BfssConfiguration& conf1, const BfssConfiguration& conf2){
Real ret(0.0);
for(size_t n(0);n<numcolor;n++){
  ret+=conf1.phases()(n)*conf2.phases()(n);
}
return ret;
}



inline Real bfssBosonicEnergy(const BfssConfiguration& conf) {
  Real coeff(-3.0/4.0);
  Real commPart(0.0);
  BfssConfiguration::MatrixType temp, comm;
  #pragma omp parallel for
    for(size_t site=0; site<conf.size(); site++){
      for(size_t e(conf.elementnum); e--;){
        for (size_t l(e); l--;) {
          temp.noalias()=conf(site, e) * conf(site, l);
          comm.noalias()=temp-temp.adjoint();
          commPart += real((comm * comm).trace());
        }
      } 
    }
  return coeff*commPart/(Real(numcolor)*conf.size());
}


inline Real correlatorX2(const BfssConfiguration& conf, int dT){
  int tPlusDt(0);
  Real corr(0.0);
  BfssConfiguration::MatrixType x2T, x2TPlusDT;
  Complex tempX2T(0.0);
  Complex tempX2TPlusDT(0.0);
  // #pragma omp parallel for
  for (size_t site=0; site<conf.size(); site++){
    tempX2T = Complex(0.0, 0.0);
    tempX2TPlusDT = Complex(0.0, 0.0);
    tPlusDt = site + dT;
    if (tPlusDt>=conf.size()){
      tPlusDt -= conf.size();
    }
    for (size_t e(conf.elementnum); e--;){
        x2T.noalias() = (conf(site, e) * conf(site, e));
        x2TPlusDT.noalias() = (conf(tPlusDt, e) * conf(tPlusDt, e));
        tempX2T+=x2T.trace();
        tempX2TPlusDT+=x2TPlusDT.trace();
    }
    corr += real(tempX2T * tempX2TPlusDT);
  }
  return corr/(conf.size());
}

inline Real correlatorX2_Corrected(const BfssConfiguration& conf, int dT){
  int tPlusDt(0);
  Real corr(0.0);
  Complex tempX2T(0.0);
  Complex tempX2TPlusDT(0.0);
  // #pragma omp parallel for
  for (size_t site=0; site<conf.size(); site++){
    tempX2T = Complex(0.0, 0.0);
    tempX2TPlusDT = Complex(0.0, 0.0);
    tPlusDt = site + dT;
    if (tPlusDt>=conf.size()){
      tPlusDt -= conf.size();
    }
    for (size_t e(conf.elementnum); e--;){
        tempX2T= (conf(site, e) * conf(site, e)).trace();
        tempX2TPlusDT = (conf(tPlusDt, e) * conf(tPlusDt, e)).trace();
        corr += real(tempX2T * tempX2TPlusDT);
    }
  }
  return corr/(conf.size());
}

inline Real corrT0(const BfssConfiguration& conf, int dT){
  return real( (conf(0, 1)*conf(dT, 1)).trace() );
}


inline Real trX2(const BfssConfiguration& conf){
  Real val(0.0);
  for (size_t site=0; site<conf.size(); site++){
    for (size_t e(conf.elementnum); e--;){
        val += real( (conf(site, e) * conf(site, e)).trace());
    }
  }
  return val/(conf.size());
}

inline Real test_corr(const BfssConfiguration& conf, int dT){
  Real corr(0.0);
  int tPlusDt(0);
  for (size_t t=0; t<conf.size(); t++){
    tPlusDt = t + dT;
    if (tPlusDt>=conf.size()){
      tPlusDt -= conf.size();
    }
    corr += real( ( ( conf(t, 0)*conf(t, 0) ).trace() )*( ( conf(tPlusDt, 0)*conf(tPlusDt, 0) ).trace() ) );
  }
  return corr/(conf.size());
}

inline Real singleMixedCorrelator(const BfssConfiguration& conf, int dT){
  Real corr(0.0);
  int tPlusDt(0);
  for (size_t t=0; t<conf.size(); t++){
    tPlusDt = t + dT;
    if (tPlusDt>=conf.size()){
      tPlusDt -= conf.size();
    }
    corr += real( ( ( conf(t, 0)*conf(t, 1) ).trace() )*( ( conf(tPlusDt, 0)*conf(tPlusDt, 1) ).trace() ) );
  }
  return corr/(conf.size());
}


#endif /* BFSSCONFIG_H_ */
