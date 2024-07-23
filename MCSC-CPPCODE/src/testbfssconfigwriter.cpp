/*
 * testbfssconfigwriter.cpp
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#define NUMCOLORS 8
#define SCALARFIELDS 9

//#include "bfssconfigwriter.h"
#include "bfsssystemstate.h"
#include "bfssconfigwriter.h"
#include "bfssbosonicaction.h"
#include "bfsshmcupdater.h"

#include <iostream>

void fillConfiguration(BfssConfiguration& inout){
  for (size_t jj = 0; jj < numcolor; ++jj) {
    inout.phases()(jj)=Real(jj)/Real(numcolor);
  }
for (size_t e = 0; e < inout.elementnum; ++e) {
for (size_t tt = 0; tt < inout.size(); ++tt)
  for (size_t ii = 0; ii < numcolor; ++ii)
    for (size_t jj = 0; jj < numcolor; ++jj) {
      Real val(Real(tt)+0.1*Real(e)+0.01*Real(ii+jj)/Real(numcolor*numcolor));
      inout(tt, e)(ii, jj)=val;
    }
}
}

Real testConfiguration(const BfssConfiguration& inout){
  Real err(0.0);
  for (size_t jj = 0; jj < numcolor; ++jj) {
    const Real diff(fabs(inout.phases()(jj)-Real(jj)/Real(numcolor)));
    if(diff>err){
      err=diff;
    }
  }
for (size_t e = 0; e < inout.elementnum; ++e) {
for (size_t tt = 0; tt < inout.size(); ++tt)
  for (size_t ii = 0; ii < numcolor; ++ii)
    for (size_t jj = 0; jj < numcolor; ++jj) {
      Real val(Real(tt)+0.1*Real(e)+0.01*Real(ii+jj)/Real(numcolor*numcolor));
      Real diff(fabs(real(inout(tt, e)(ii, jj))-val));
      if(diff>err){
        err=diff;
      }
    }
}
return err;
}



int main(){
  std::cout<<"Starting test of basic configuration for BFSS model"<<std::endl;
  BfssSystemState state(100);
  fillConfiguration(state.config());
  ConfigWriter writer;
  writer.initialize(state);
  std::vector<Real> data({1.0,2.0,3.0});
  std::list<std::string> names({"one","two","three"});
  writer.writeConfigurationToOctaveFormat(state,data,names);
  BfssSystemState state2(100);
  writer.readConfigurationFromOctaveFormat(state2);
  std::cout<<"Error "<<testConfiguration(state2.config())<<std::endl;
  state2.parameters().confignumber=1;
  writer.writeConfigurationToOctaveFormat(state2,data,names);
  return 0;
}



