/*
 * bechmark.cpp
 *
 *  Created on: 06.04.2023
 *      Author: Georg Bergner
 */
//#define DYNAMIC_ALLOCATION
#define NUMCOLORS 4
#define SCALARFIELDS 9

#include "basicdef.h"
#include "bfsssystemstate.h"
#include "diracoperator.h"

#include <chrono>
#include <iostream>
#include <omp.h>

const unsigned long runs(1000);

int main(){
  omp_set_num_threads(1);
  std::cout<<"Number of colors "<<numcolor<<" number of scalars "<<SCALARFIELDS<<std::endl;
  std::cout << "Executing with OMP " << omp_get_max_threads()
      << " max threads." << std::endl;
  BfssSystemState state(8);
  state.parameters().temperature=1.0;
  state.parameters().diracParameter.fermmass=0.0;
  state.parameters().diracParameter.mu=0.0;
  state.parameters().diracParameter.improved=true;
  state.parameters().diracParameter.setMasanorisConventions(state.parameters().latticeSpacing());
  std::cout<<"lattice size "<<state.config().size()<<" mu "<<state.parameters().diracParameter.mu<<" improved "<<state.parameters().diracParameter.improved<<std::endl;
  BfssDiracVector in(state.config().size());
  BfssDiracVector out(state.config().size());
  auto dirac=state.getDiracOperator();

  using std::chrono::high_resolution_clock;
   using std::chrono::duration_cast;
   using std::chrono::duration;
   using std::chrono::milliseconds;
   using std::chrono::seconds;

   auto t1 = high_resolution_clock::now();
  for(unsigned long i(0);i<runs;i++){
     dirac.multiply(in,out);
  }
  auto t2 = high_resolution_clock::now();
  std::cout<<"S "<< duration_cast<seconds>(t2 - t1).count()<<std::endl;
  std::cout<<"MS "<< duration_cast<milliseconds>(t2 - t1).count()<<std::endl;
  std::cout<<"MS/Mult "<< double(duration_cast<milliseconds>(t2 - t1).count())/double(runs)<<std::endl;
}


