/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

//#define DYNAMIC_ALLOCATION
#define NUMCOLORS 3
#define SCALARFIELDS 9

#include "simulationrunner.h"
#include "bfsssystemstate.h"
#include "bfsshmcupdater.h"
#include "bfssmeasurements.h"
#include "configparser.h"

#include <iostream>

int main(int argc, char** argv) {
  typedef SimulationRunner SRunner;
  SRunner runner;
  {
    Parser::ConfigParser reader;
    reader.read("config.cfg");
    int numthreads=reader.get_value_noerr("Machine","Threads",int(0));
    if(numthreads>0){
      omp_set_num_threads(numthreads);
    }
    LOG_STREAM<<"Number of colors "<<numcolor<<" number of scalars "<<SCALARFIELDS<<std::endl;
    LOG_STREAM << "Executing with OMP " << omp_get_max_threads()
        << " max threads." << std::endl;
    runner.initialize(reader);
  }
  runner.run();
  return 0;
}
