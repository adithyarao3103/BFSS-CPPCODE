/*
 * testbosonicaction.cpp
 *
 *  Created on: 24.03.2023
 *      Author: Georg Bergner
 */

#define NUMCOLORS 4
#define SCALARFIELDS 9

#include "bfsssystemstate.h"
#include "bfsshmcupdater.h"

#include <iostream>

void runTest(){
  std::cout<<"Printing result of some sample configuration to compare."<<std::endl;
  BfssSystemState state(8);
  state.parameters().temperature=0.5;
  state.parameters().bosmass=0.0;
  state.parameters().bosonmasssq=0.0;
  state.parameters().mu=4.0;
  state.parameters().improvedaction=true;
  state.config().setToZero();
  state.config().setPhaseToIdentity();
  for(size_t i(0);i<numcolor;i++){
    state.config().phases()(i)=0.4*Real(i+1)*Real(i+1)-1.5;
  }
  state.config()(0,0)(0,1)=1.0;
  state.config()(0,0)(1,0)=1.0;
  state.config()(0,1)(0,0)=1.0;
  state.config()(0,1)(1,1)=-1.0;
  state.config()(0,2)(0,1)=Complex(0.0,-1.0);
  state.config()(0,2)(1,0)=Complex(0.0,1.0);
  state.config()(1,0)(0,1)=1.0;
  state.config()(1,0)(1,0)=1.0;
  state.config()(2,0)(0,1)=1.0;
  state.config()(2,0)(1,0)=1.0;
  std::cout<<"TrX2 "<<sumTraceX2Norm(state.config());
  std::cout<<" PL "<<polyakovLoopAlpha(state.config())<<std::endl;
  std::cout<<" kin prefact "<<0.5 * static_cast<Real>(numcolor) *state.parameters().kinPrefact/state.parameters().latticeSpacing()<<std::endl;
  BfssBosonicAction bosonaction;
  bosonaction.update(state);
  bosonaction.logActionParts(std::cout);
  std::cout<<std::endl;
}

void testRandomFields(BfssSystemState& state){
  std::cout<<"Basic random number generation test."<<std::endl;
  RandomNumberGenerator rand;
  setGaussianRandomX(state.config(),rand);
  Real averalpha=0.0;
  for(int num=0;num<100;num++){
    setGaussianRandomAlpha(state.config(),rand);
    averalpha+=sumAlpha2(state.config())/Real(numcolor);
  }
  averalpha/=100;
  std::cout<<"average X2 "<<sumTraceX2(state.config())/Real(numcolor*numcolor*BfssConfiguration::elementnum*state.config().size())<<" should be 0.5 "<<std::endl;
  std::cout<<"average Alpha2 "<<averalpha<<" should be 0.5 "<<std::endl;
}

int main(){
  runTest();
  BfssSystemState state(20);
  testRandomFields(state);
}



