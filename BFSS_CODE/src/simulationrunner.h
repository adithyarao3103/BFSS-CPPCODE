/*
 * simulationrunner.h
 *
 *  Created on: 24.03.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_SIMULATIONRUNNER_H_
#define SRC_BFSS_SIMULATIONRUNNER_H_

#include "basicdef.h"
#include "configparser.h"
#include "bfsshmcupdater.h"
#include "bfsssystemstate.h"
#include "bfssconfigwriter.h"
#include "bfssmeasurements.h"

#include <iostream>
#include <ctime>

class Time{
private:
  time_t startt_={};
public:
  Time():startt_(){
    startt_=time(0);
  }

  std::string current()const {
      time_t     now = time(0);
      struct tm  tstruct;
      char       buf[80];
      tstruct = *localtime(&now);
      double seconds=difftime(now,startt_);
      strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
      std::ostringstream os;
      os<<"("<<int(seconds)<<")";
      return buf+os.str();
  }


};




struct SimulationParameter{
    unsigned long updateSteps=0;
    unsigned long warmupSteps=0;
    unsigned long measureEach=1;
    unsigned long startConfig=0;
    bool readstartconfig=false;
    unsigned int writeEach=0;
    bool continueWithLast=false;
};



inline void readSimulationParameter(Parser::ConfigParser& parser,SimulationParameter& param){
  param.updateSteps=parser.get_value_noerr("Simulation","UpdateSteps",param.updateSteps);
  param.warmupSteps=parser.get_value_noerr("Simulation","WarmUpSteps",param.warmupSteps);
  param.measureEach=parser.get_value_noerr("Simulation","MeasureEach",param.measureEach);
  param.startConfig=parser.get_value_noerr("Simulation","StartConfig",param.startConfig);
  param.readstartconfig=parser.get_value_noerr("Simulation","ReadStart",param.readstartconfig);
  param.writeEach=parser.get_value_noerr("Simulation","WriteEach",param.writeEach);
  param.continueWithLast=parser.get_value_noerr("Simulation","Continue",param.continueWithLast);
  LOG_STREAM<<"Simulation:UpdateSteps:"<<param.updateSteps<<":WarmUpSteps:"
      <<param.warmupSteps<<":MeasureEach:"<<param.measureEach<<":StartConfig:"
      <<param.measureEach<<":ReadStart:"<<param.readstartconfig
      <<":WriteEach:"<<param.writeEach<<":Continue:"<<param.continueWithLast<<std::endl;
}


class SimulationRunner {
public:
    using Measurement=BfssMeasurements;
    using Updater=BfssUpdaterHmc;
    using StateFactory=BfssSystemStateFactory;
    using State=StateFactory::State;
    using Reader=ConfigWriter;

protected:
  State* state_=nullptr;
  Updater updater_={};
  Measurement measurement_={};
  Reader readwriter_={};
  SimulationParameter parameter_={};
  Time time_={};

public:
  SimulationRunner(const SimulationRunner&)=delete;
  SimulationRunner& operator=(const SimulationRunner&)=delete;

  SimulationRunner()
    {
  }

  void initialize(Parser::ConfigParser& reader) {
    printTiming_(" Initializing Simulation");
    readSimulationParameter(reader,parameter_);
    StateFactory fact;
    state_ = fact.getState(reader);
    updater_.initialize(*state_, reader);
    measurement_.initialize(*state_, reader);
    measurement_.printStart();
    printTiming_(" Finished Initializing Simulation");
  }

  virtual ~SimulationRunner() {
    delete state_;
  }

  virtual void run() {
    printTiming_(" Starting Initialize Simulation");
    initializeState_(*state_);
    updater_.runChecks(*state_);
    printTiming_(" Starting Warmup Simulation");
    runWarmUp_(*state_);
    updater_.reset();
    printTiming_(" Starting Update Simulation");
    runUpdate_(*state_);
    printTiming_(" End Update Simulation");
    finalize_(*state_);
  }

private:

   void initializeState_(State& state){
     state.configNum() = parameter_.startConfig;
     updater_.initializeFields(state);
     bool found=false;
     if(parameter_.continueWithLast){
       std::ifstream conflog("writtenconfigs.log");
       unsigned long r;
       while(conflog){
         conflog>>r;
         found=true;
       }
       if(found){
        LOG_STREAM<<"Found last written config "<<r<<std::endl;
        state.configNum()=r;
       }
     }
     if(parameter_.readstartconfig||parameter_.continueWithLast){
       readwriter_.initialize(state);
       if(found||parameter_.readstartconfig){
        readwriter_.readConfig(state);
       }
     }else{
       if(parameter_.writeEach>0){
         readwriter_.initialize(state);
       }
     }
   }

   void printTiming_(const std::string& str)const{
     LOG_STREAM<<time_.current()<< str<<std::endl;
   }



  void runWarmUp_(State& state){
    for (unsigned int i(0); i < parameter_.warmupSteps; i++ ) {
      updater_.update(state);
      state.configNum()++;

    }
  }

   void runUpdate_(State& state) {
    unsigned long meas(0);
    unsigned long cwrite(0);
    for (size_t up(1); up <= parameter_.updateSteps; up++) {
      updater_.update(state);
      if (parameter_.measureEach > 0 && up % parameter_.measureEach == 0) {
        measurement_.measure(state);
        meas++;
      }
      if(parameter_.writeEach>0 && up%parameter_.writeEach==0){
        readwriter_.writeConfig(state);
        std::ofstream conflog("writtenconfigs.log",std::ios_base::app);
        conflog<<state.configNum()<<std::endl;
        cwrite++;
      }
      state.configNum()++;
    }
    LOG_STREAM<<"Number of measurements "<<meas<<" number of written configs "<<cwrite<<std::endl;
    }

  void finalize_(State& state) {
    measurement_.printFinal();
    updater_.printData(std::cout, state);
  }



};



#endif /* SRC_BFSS_SIMULATIONRUNNER_H_ */
