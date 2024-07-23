/*
 * testdiracoperator.cpp
 *
 *  Created on: 27.03.2023
 *      Author: Georg Bergner
 */

#define NUMCOLORS 4
#define SCALARFIELDS 9

#include "diracoperator.h"
#include "bfsssystemstate.h"
#include "cgmsolver.h"

template<class T>
void print10DGammaStructure(){
    for(size_t al1=0;al1<16;al1++){
      for(size_t al2=0;al2<16;al2++){
        std::cout<<(al2==T::index(al1) ? T::value(al1):0)<<" ";
      }
      std::cout<<"\n";
    }
    std::cout<<"\n";
}

template<unsigned int mu, bool minus=false>
void testSquare10DGammaStructure(){
   using T=Gamma10D<mu>;
    for(size_t al1=0;al1<16;al1++){
      const size_t ind=T::index(T::index(al1));
      const Complex val=T::value(T::index(al1))*T::value(al1);
      if(ind!=al1){
        std::cout<<"Error sq with index: "<<ind<<"!="<<T::index(T::index(al1))<<std::endl;
      }
      if(std::abs(val-(minus ? -Real(1.0):Real(1.0)))>1e-3){
        std::cout<<"Error sq with value: "<<val<<std::endl;
      }
    }
}

template<unsigned int mu,unsigned int nu>
void testAnticommutator10DGammaStructure(){
  using T1=Gamma10D<mu>;
  using T2=Gamma10D<nu>;
    for(size_t al1=0;al1<16;al1++){
      const size_t ind1=T1::index(T2::index(al1));
      const size_t ind2=T2::index(T1::index(al1));
      const Complex val1=T1::value(T2::index(al1))*T2::value(al1);
      const Complex val2=T2::value(T1::index(al1))*T1::value(al1);
      if(ind1!=ind2){
        std::cout<<"Error aComm "<<mu<<" "<<nu<<" with index: "<<ind1<<"!="<<ind2<<std::endl;
      }
      if(std::abs(val1+val2)>1e-3){
        std::cout<<"Error aComm "<<mu<<" "<<nu<<" with value: "<<val1<<" "<<val2<<std::endl;
      }
    }
}

void testGamma(){
  //  print10DGammaStructure<Gamma10D<4>>();
    testSquare10DGammaStructure<0>();
    testSquare10DGammaStructure<1>();
    testSquare10DGammaStructure<2>();
    testSquare10DGammaStructure<3>();
    testSquare10DGammaStructure<4>();
    testSquare10DGammaStructure<5>();
    testSquare10DGammaStructure<6>();
    testSquare10DGammaStructure<7>();
    testSquare10DGammaStructure<8,true>();
    testSquare10DGammaStructure<9>();
    //  testAnticommutator10DGammaStructure<0,8>(); This should of course not work.
    testAnticommutator10DGammaStructure<0,1>();
    testAnticommutator10DGammaStructure<0,2>();
    testAnticommutator10DGammaStructure<0,3>();
    testAnticommutator10DGammaStructure<0,4>();
    testAnticommutator10DGammaStructure<0,5>();
    testAnticommutator10DGammaStructure<0,6>();
    testAnticommutator10DGammaStructure<0,7>();
    testAnticommutator10DGammaStructure<0,9>();
    testAnticommutator10DGammaStructure<1,2>();
    testAnticommutator10DGammaStructure<1,3>();
    testAnticommutator10DGammaStructure<1,4>();
    testAnticommutator10DGammaStructure<1,5>();
    testAnticommutator10DGammaStructure<1,6>();
    testAnticommutator10DGammaStructure<1,7>();
    testAnticommutator10DGammaStructure<1,9>();
    testAnticommutator10DGammaStructure<2,3>();
    testAnticommutator10DGammaStructure<2,4>();
    testAnticommutator10DGammaStructure<2,5>();
    testAnticommutator10DGammaStructure<2,6>();
    testAnticommutator10DGammaStructure<2,7>();
    testAnticommutator10DGammaStructure<2,9>();
    testAnticommutator10DGammaStructure<3,4>();
    testAnticommutator10DGammaStructure<3,5>();
    testAnticommutator10DGammaStructure<3,6>();
    testAnticommutator10DGammaStructure<3,7>();
    testAnticommutator10DGammaStructure<3,9>();
    testAnticommutator10DGammaStructure<4,5>();
    testAnticommutator10DGammaStructure<4,6>();
    testAnticommutator10DGammaStructure<4,7>();
    testAnticommutator10DGammaStructure<4,9>();
    testAnticommutator10DGammaStructure<5,6>();
    testAnticommutator10DGammaStructure<5,7>();
    testAnticommutator10DGammaStructure<5,9>();
    testAnticommutator10DGammaStructure<6,7>();
    testAnticommutator10DGammaStructure<6,9>();
    testAnticommutator10DGammaStructure<7,9>();
}

class FakeOperator {
  public:
    using Vector=BfssDiracVector;
  void multiply(const Vector& in,Vector& out){
#pragma omp parallel for
     for (size_t site = 0; site <in.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         out(site,e)=Complex(0.0,0.5)*in(site,e);
       }
     }
#pragma omp parallel for
     for (size_t site = 0; site <in.size(); site++) {
       const size_t up = in.up(site);
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         out(up,e)+=Complex(0.0,1.5)*in(site,e);
       }
     }
     out(3,0)+=0.1*in(4,0);
     out(0,2)+=0.1*in(0,2);
     out(1,0)+=0.1*in(1,0);
  }
};

class FakeOperatorDag {
  public:
    using Vector=BfssDiracVector;
  void multiply(const Vector& in,Vector& out){
#pragma omp parallel for
     for (size_t site = 0; site <in.size(); site++) {
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         out(site,e)=Complex(0.0,-0.5)*in(site,e);
       }
     }
#pragma omp parallel for
     for (size_t site = 0; site <in.size(); site++) {
       const size_t up = in.up(site);
       for (size_t e(BfssDiracVector::elementnum); e--;) {
         out(site,e)+=-Complex(0.0,1.5)*in(up,e);
       }
     }
     out(4,0)+=0.1*in(3,0);
     out(0,2)+=0.1*in(0,2);
     out(1,0)+=0.1*in(1,0);
  }
};

void testOutputComparison(){
  BfssSystemState state(8);
  state.parameters().temperature=0.5;
  state.parameters().diracParameter.mu=4.0;
  state.parameters().diracParameter.improved=true;
  state.parameters().diracParameter.fermionbc=-1.0;
  state.parameters().diracParameter.setMasanorisConventions(state.parameters().latticeSpacing());
  auto dirac=state.getDiracOperator();
  BfssDiracVector in(8);
  BfssDiracVector out(8);
  state.config().setToZero();
  in.setToZero();
  out.setToZero();
  for(size_t i(0);i<numcolor;i++){
   state.config().phases()(i)=0.4*Real(i+1)*Real(i+1)-1.5;
  }
//  state.config()(0,1)(0,1)=0.1;
 // state.config()(0,1)(1,0)=0.1;
  state.config()(0,8)(0,1)=1.1;
  state.config()(0,8)(1,0)=1.1;
  for(size_t i(0);i<in.elementnum;i++){
    in(0,i)(0,0)=1.0+i;
    in(1,i)(1,0)=Complex(0.0,2.1);
    in(0,i)(0,0)=1.0+i+Complex(0.0,0.5);
  }
  in(0,8)(0,0)=1.0;
  dirac.multiply(in,out);
  std::cout<<"Prefact kin "<<state.parameters().diracParameter.prefactorKin<<" Prefact comm "<<state.parameters().diracParameter.prefactComm<<std::endl;
  std::cout<<"In vector "<<normVectorSq(in)<<" out vector "<<normVectorSq(out)<<" polyakov "<<polyakovLoopAlpha(state.config())<<" trx2 "<<sumTraceX2Norm(state.config())<<std::endl;
  std::cout<<"site=0 --- "<<out(0,0)(0,0)<<" "<<out(0,8)(0,0)<<" "<<out(1,0)(0,0)<<" "<<out(1,8)(0,0)<<std::endl;
  std::cout<<"site=L-1 --- "<<out(out.size()-1,0)(0,0)<<" "<<out(out.size()-1,8)(0,0)<<std::endl;
}

void testDDagger(){
  BfssSystemState state(8);
  state.parameters().temperature=0.5;
  state.parameters().diracParameter.mu=4.0;
  state.parameters().diracParameter.improved=true;
  state.parameters().diracParameter.fermionbc=-1.0;
  state.parameters().diracParameter.setMasanorisConventions(state.parameters().latticeSpacing());
  BfssDiracVector in(8);
  BfssDiracVector out(8);
  BfssDiracVector outdag(8),in2(8);
  RandomNumberGenerator rnd;
  setGaussianRandomVector(in,rnd);
  setGaussianRandomVector(in2,rnd);
  setGaussianRandomX(state.config(),rnd);
  auto dirac=state.getDiracOperator();
  dirac.multiply(in,out);
  dirac.setDag(true);
  dirac.multiply(in,outdag);
  std::cout<<"<v,D v> "<<scalarProd(in,out)<<"<D^dagv,v> "<<scalarProd(outdag,in)<<std::endl;
  dirac.multiply(in2,outdag);
  std::cout<<"<v,D w> "<<scalarProd(in2,out)<<"<D^dagv,w> "<<scalarProd(outdag,in)<<std::endl;
}

void testInverters(){
  BfssSystemState state(8);
  BfssDiracVector in(8);
  BfssDiracVector out(8);
  BfssDiracVector outdag(8),in2(8);
  RandomNumberGenerator rnd;
  setGaussianRandomVector(in,rnd);
  setGaussianRandomVector(in2,rnd);
  setGaussianRandomX(state.config(),rnd);
  state.parameters().diracParameter.prefactorKin=1.0;
  state.parameters().diracParameter.prefactComm=1.0;
  state.parameters().diracParameter.prefactMu=1.0;
  state.parameters().diracParameter.improved=true;
  state.parameters().diracParameter.mu=10.0;
  state.parameters().diracParameter.setMasanorisConventions(state.parameters().latticeSpacing());
  auto dirac=state.getDiracOperator();
  dirac.setDag(true);
  CgnrSolver<BfssDiracOperator,BfssDiracOperator> solver(state.getDiracOperator(),dirac);
  CgSqSolver<BfssDiracOperator,BfssDiracOperator> solver2(state.getDiracOperator(),dirac);
  solver.setMaxIterations(100000);
  solver.setResidGoal(1e-5);
  out.setToZero();
  projectTracelessVector(in);
  solver.multiply(in,out);
  std::cout<<"Checking trace "<<sumAbsTraceVector(in)<<" "<<sumAbsTraceVector(out)<<std::endl;
  dirac.setDag(false);
  dirac.multiply(out,in2);
  rescaleAdd(in2,-1.0,in);
  std::cout<<"Check resid "<<sqrt(normVectorSq(in2)/normVectorSq(in))<<std::endl;
  solver2.setMaxIterations(1000);
  solver2.setResidGoal(1e-5);
  out.setToZero();
  projectTracelessVector(in);
  solver2.multiply(in,out);
  dirac.setDag(true);
  dirac.multiply(out,outdag);
  dirac.setDag(false);
  dirac.multiply(outdag,in2);
  rescaleAdd(in2,-1.0,in);
  std::cout<<"Checking trace "<<sumAbsTraceVector(in)<<" "<<sumAbsTraceVector(out)<<std::endl;
  std::cout<<"Check resid "<<sqrt(normVectorSq(in2)/normVectorSq(in))<<std::endl;
  dirac.setDag(true);
  std::vector<Real> shifts({0.1,0.2,0.3,0.4});
  CgSqMSolver<BfssDiracOperator,BfssDiracOperator> solverm(state.getDiracOperator(),dirac,shifts);
  std::vector<BfssDiracVector> outv(shifts.size());
  solverm.multiply(in,outv);
  for(size_t i(0);i<outv.size();i++){
    dirac.setDag(false);
   dirac.multiply(outv[i],outdag);
   dirac.setDag(true);
   dirac.multiply(outdag,in2);
   addRescaled(in2,outv[i],shifts[i]);
   rescaleAdd(in2,-1.0,in);
   std::cout<<"Check resid shift "<<i<<" "<<sqrt(normVectorSq(in2)/normVectorSq(in))<<std::endl;
  }

}

int main(){
  testGamma();
  testOutputComparison();
  testDDagger();
  testInverters();
}
