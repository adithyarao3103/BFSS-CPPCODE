#ifndef FERMION_OBSERVABLES_H_
#define FERMION_OBSERVABLES_H_

#include "basicdef.h"
#include "cgmsolver.h"
#include "diracoperator.h"
#include "bfsssystemstate.h"
#include "bfssconfig.h"
#include "RandomNumberGenerator.h"

using Vector = BfssDiracVector;
using DiracOperator = BfssDiracOperator;

template<int n>
struct loopCommutator{
    static void exec(const size_t site, Vector &in, Vector &out, const BfssConfiguration &conf){
        for (size_t elem(0); elem < Vector::elementnum; elem++){
            out(site, elem).noalias() += Gamma10D<n>::value(elem) * ( conf(site, n)*in(site, Gamma10D<n>::index(elem)) - in(site, Gamma10D<n>::index(elem))*conf(site, n) );
        }
        loopCommutator<n-1>::exec(site, in, out, conf);
    }
};

template<>
struct loopCommutator<-1>{
    static void exec(const size_t site, Vector &in, Vector &out, const BfssConfiguration &conf){

    }
};

/*
From Masanori's manual,  Fermionic energy = < 3/2 Tr(D^-1 K) / N\beta > where K = \gamma^M [X_M , . ]
This can be calculated as < 3/2 phi^dagger D^-1 K phi / (phi^dagger phi) N \beta >
*/

inline Complex bfssFermionicEnergy(const BfssConfiguration &conf, int numNoise, const DiracOperator &op, const DiracOperator &opdag){
    Complex energy(0.0, 0.0);
    RandomNumberGenerator rnd;
    CgnrSolver<DiracOperator, DiracOperator> solver(op, opdag);
    // Neccessary Vectors
    Vector phi(conf.size()), KPhi(conf.size()), DInvKPhi(conf.size());
    
    for (int i = 0; i< numNoise; i++){
        setGaussianRandomVector(phi, rnd);
        // Act by K
        for(size_t site = 0; site< conf.size(); site++){
            loopCommutator<BfssConfiguration::elementnum - 1>::exec(site, phi, KPhi, conf);
        }
        // Act by DInv
        solver.multiply(KPhi, DInvKPhi);
        // Find inner product with original noise vector
        Complex inner(0.0, 0.0);
        for(size_t site = 0; site< conf.size(); site++){
            for(size_t elem(Vector::elementnum);elem--;){
                inner += (phi(site, elem).adjoint() * DInvKPhi(site, elem)).trace();
            }
        }
        // To remove the historical i in front of the Dirac Operator
        inner = inner*Complex(0, 1);
        energy += inner;
    }
    //Normalisation: There should be an extra 2 and numcolor in the denominator due to the wrong normalisation of the Dirac Operator. 
    //The extra two in the denominator is cancelled by the one in numerator which arises due to the extra division by 2 in the setGaussianRandomVector function.
    Real prefact(3.0/(2.0*numcolor*numcolor*conf.size()));
    return Complex(energy.real()*prefact/numNoise, energy.imag()*prefact/numNoise);
}


inline Complex traceOfIdentity(const BfssConfiguration &conf, int numNoise){
    Complex trace(0.0, 0.0);
    RandomNumberGenerator rnd;
    const size_t latticeSize = conf.size();
    // Neccessary Vectors
    Vector phi(latticeSize);
    
    for (int i = 0; i< numNoise; i++){
        setGaussianRandomVector(phi, rnd);
        for(size_t site = 0; site< conf.size(); site++){
            for(size_t elem(Vector::elementnum);elem--;){
                trace += (phi(site, elem).adjoint() * phi(site, elem)).trace();
            }
        }
    }
    //Normalisation 
    return Complex(trace.real()*2/numNoise, trace.imag()*2/numNoise);
}



#endif
