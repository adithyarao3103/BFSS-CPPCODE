/*
 *
 *  Created on: 24.03.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_BFSSMEASUREMENTS_H_
#define SRC_BFSS_BFSSMEASUREMENTS_H_
#include "bfsssystemstate.h"
#include "configparser.h"
#include "fermionmeasurements.h"

#include <fstream>

  class BfssMeasurements{
  private:  
   std::ofstream outputfile_={};
  public:
    using SystemState=BfssSystemState;

    void initialize(SystemState& state,
        Parser::ConfigParser& parser) {
       outputfile_.open("measurements.dat",std::ios_base::app);
    }

    void measure(SystemState& state){
      outputfile_<<state.configNum()<<" ";
      Complex ferm(bfssFermionicEnergy(state.config(), 20, state.getDiracOperator(), state.getDiracOperatorDagger()));
      // Complex ferm(0,0);
      Real bos(bfssBosonicEnergy(state.config()));
      outputfile_<<ferm.real()<<" "<<ferm.imag()<<" ";
      outputfile_<<bos<<" ";
      outputfile_<<ferm.real() + bos<<" ";
      outputfile_<<sumTraceX2Norm(state.config())<<" ";
      outputfile_<<traceOfIdentity(state.config(), 20);
      outputfile_<<std::endl;
    }

    void printStart(){
      outputfile_<<"# NUMCOLORS "<<Real(numcolor)<<" ";
      outputfile_<<"Total Energy Functional";
      outputfile_<<std::endl;
      outputfile_<<"#1:ConfigNumber ";
      outputfile_<<"2:Fermionic Energy (Real) ";
      outputfile_<<"3:Fermionic Energy (Imag) ";
      outputfile_<<"4:Bosonic Energy ";
      outputfile_<<"5:Total Energy ";
      outputfile_<<"6:X^2 ";
      outputfile_<<"7:Tr(I)";
      outputfile_<<std::endl;
    }

    void printFinal(){

    }

}; 

#endif /* SRC_BFSS_BFSSMEASUREMENTS_H_ */
