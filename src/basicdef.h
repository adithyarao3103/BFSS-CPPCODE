/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef BASICDEF_H
#define BASICDEF_H

#define EIGEN_NO_DEBUG
// Math
#include <cstring>
#include <complex>
#include <vector>
#include <cmath>

// Eigen
#include <Core>
#include <Dense>
#include <Sparse>

using Real=double;

#include <boost/math/constants/constants.hpp>
const Real pi=boost::math::constants::pi<Real >();// since we can not assume C++20 for now

#ifndef NUMCOLORS
 #define NUMCOLORS 4
 #warning "NUMCOLORS not set using default 4"
#endif

#ifndef SCALARFIELDS
  #define SCALARFIELDS 9
  #warning "SCALARFIELDS not set using default 9"
#endif

#ifndef BFSSSPINSIZE
 #define BFSSSPINSIZE 16
#endif

#if NUMCOLORS==2
#define SU2
#endif



#if NUMCOLORS>12
#ifndef DYNAMIC_ALLOCATION
#define DYNAMIC_ALLOCATION
#endif
#endif

const unsigned int numcolor = NUMCOLORS;
const unsigned int gaugerank = numcolor*numcolor-1;


using Complex=std::complex<Real>;
using GaugeMatrix=Eigen::Matrix<Complex, numcolor, numcolor>;
using GaugeVector=Eigen::Matrix<Complex, numcolor,1>;
using DiagonalGaugeMatrix=Eigen::DiagonalMatrix<Complex,numcolor>;
using RealGaugeVector=Eigen::Matrix<Real, numcolor,1>;
using DiagonalMatrix=Eigen::DiagonalMatrix<Complex, numcolor, numcolor>;
using DynamicGaugeMatrix=Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic >;
using GaugeMatrixMap=Eigen::Map<DynamicGaugeMatrix >;
using GaugeMatrixMapConst=Eigen::Map<const DynamicGaugeMatrix >;
using std::conj;


// Helper functions not provided by Eigen library:
auto maxIndex(const RealGaugeVector& v){
   // does not catch case of v.size()==0
    size_t ind(0);
    Real mval=v[0];
    for(size_t i(0);i<v.size();i++){
      if(v[i]>mval){
        ind=i;
        mval=v[i];
      }
    }
    return std::pair<Real,size_t>(mval,ind);
}

auto minIndex(const RealGaugeVector& v){
   // does not catch case of v.size()==0
    size_t ind(0);
    Real mval=v[0];
    for(size_t i(0);i<v.size();i++){
      if(v[i]<mval){
        ind=i;
        mval=v[i];
      }
    }
    return std::pair<Real,size_t>(mval,ind);
}

// default matrix generation
template<class T>
inline T generateMatrix(){
  return GaugeMatrix();
}


template<>
inline DynamicGaugeMatrix generateMatrix<DynamicGaugeMatrix>(){
  return DynamicGaugeMatrix(numcolor,numcolor);
}

#define LOG_STREAM std::cout
#define LOG_MESSAGE(X) std::cout<< X <<std::endl;

#endif
