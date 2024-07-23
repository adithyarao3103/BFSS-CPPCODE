/*
 * testrationalapproximation.cpp
 *
 *  Created on: 05.04.2023
 *      Author: Georg Bergner
 */

#define NUMCOLORS 2
#define SCALARFIELDS 9

#include "bfssdiracaction.h"
#include <iostream>
#include <sstream>
#include <fstream>


int main(int argc, char** argv){
  if(argc<5){
    std::cout<<"Use arguments: filename, minev, maxev, steps"<<std::endl;
    return 0;
  }
  std::string filename(argv[1]);
  Real minev(0);
  Real maxev(100);
  unsigned long steps(10);
  {
   std::istringstream iss(argv[2]);
   iss>>minev;
  }
  {
   std::istringstream iss(argv[3]);
   iss>>maxev;
  }
  {
   std::istringstream iss(argv[4]);
   iss>>steps;
  }
  std::cout<<"Plot range: "<<minev<<" to "<<maxev<<" with steps "<<steps<<std::endl;
  const Real delta((maxev-minev)/Real(steps));
  RationalApproximation rat;
  readRationalApproximation(rat,filename,0);
  if(rat.shiftSize()!=0){
    rat.print(std::cout);
    std::ofstream of("testrational1.txt");
    Real x(minev);
    while(x<maxev){
      of<<x<<" "<<rat(x)<<"\n";
      x+=delta;
    }
  }
  readRationalApproximation(rat,filename,1);
  if(rat.shiftSize()!=0){
    rat.print(std::cout);
    Real x(minev);
    std::ofstream of("testrational2.txt");
    while(x<maxev){
      of<<x<<" "<<rat(x)<<"\n";
      x+=delta;
    }
  }
  return 0;
}
