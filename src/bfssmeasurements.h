/*
 *
 *  Created on: 24.03.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_BFSSMEASUREMENTS_H_
#define SRC_BFSS_BFSSMEASUREMENTS_H_
#include "bfsssystemstate.h"
#include "configparser.h"

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
      outputfile_<<sumTraceX2Norm(state.config())<<" ";
      outputfile_<<polyakovLoopAlpha(state.config())<<" ";
      outputfile_<<std::endl;
    }

    void printStart(){
      outputfile_<<"#1:ConfigNumber ";
      outputfile_<<"2:X2 ";
      outputfile_<<"3:PL ";
      outputfile_<<std::endl;
    }

    void printFinal(){

    }

};



#endif /* SRC_BFSS_BFSSMEASUREMENTS_H_ */
