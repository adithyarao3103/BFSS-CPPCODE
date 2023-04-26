/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_BFSSSYSTEMSTATE_H_
#define SRC_BFSS_BFSSSYSTEMSTATE_H_

#include "bfssconfig.h"
#include "diracoperator.h"
#include "configparser.h"

#include <map>



struct BfssParameterSet{
    Real temperature=1.0;
    size_t latticeSize=0;

    // Simulation technical parameters:
    // The following avoids NaN values in log.
    Real gaugefixRegulator=1e-8;
    Real alphaConstraintRegulator=1e-8;
    Real alphaConstraintRepulsion=100;
    // constraint on flat directions
    Real trX2Cut=0.0;
    Real trX2CutG=0.0;

    // to deviate from standard Bfss
    Real bosmass=0.0;
    Real bosonmasssq=0.0; // allows for negative square like in Higgs models.
    bool ungauged=false;

    // not implemented
    Real fermionbc=1.0;
    bool improvedaction=false;
    Real lambda=0.0;
    Real mu=0.0;
    Real fermmass=0.0;

    // These pre-factors are used mainly for debugging to deactivate parts of the action.
    Real commPrefact=1.0;
    Real kinPrefact=1.0;
    Real alphaConstraintPrefact=1.0;
    Real alphaFixPrefact=1.0;
    bool backroundgauge=false;

    // updated by the simulation runner.
    size_t confignumber=0;

    Real bosonMassSquare()const{
      return bosmass*bosmass+bosonmasssq;
    }

    Real latticeSpacing()const{
      return  1.0/(temperature*static_cast<Real>(latticeSize));
    }
    BfssDiracOperatorParameter diracParameter;
};

inline void readStateParameter(Parser::ConfigParser& parser, BfssParameterSet& param){
  // can not use defaults for the following parameters:
  param.temperature=parser.get_value<Real>("Model","Temperature");
  param.latticeSize=parser.get_value<size_t>("Model","LatticeSize");
  param.mu=parser.get_value_noerr("Model","Mu",Real(0.0));
  param.fermmass=parser.get_value_noerr("Model","FermMass",Real(0.0));
  param.fermionbc=parser.get_value_noerr("Model","FermionBC",param.fermionbc);
  param.bosonmasssq=parser.get_value_noerr("Model","BosonMassSq",param.bosonmasssq);
  param.improvedaction=parser.get_value_noerr("Model","Improved",param.improvedaction);
  param.alphaConstraintRepulsion=parser.get_value_noerr("Model","AlphaConstraintRepulsion",param.alphaConstraintRepulsion);
  LOG_STREAM<<"Model:Temperature:"<<param.temperature<<":LatticeSize:"<<param.latticeSize<<":Mu:"<<param.mu<<":FermionBC:"<<param.fermionbc<<":BosonMassSq:"<<param.bosonmasssq
      <<":FermMass:"<<param.fermmass<<":Improved:"<<param.improvedaction<<std::endl;
  LOG_STREAM<<":AlphaConstraintRepulsion:"<<param.alphaConstraintRepulsion<<std::endl;
}

struct StoredDataSet{
    Real bosonicAction=0.0;
    Real fermionicAction=0.0;
    Real bosonicKineticPart=0.0;
    Real bosonicCommutatorPart=0.0;
    Real bosonicBMNPart=0.0;
    Real momentumPart=0.0;
    Complex polyakovLoop=0.0;
};

class BfssSystemState {
private:
  BfssParameterSet parameters_={};
  StoredDataSet storedDataSet_={};
  BfssConfiguration configuration_;

public:

  explicit BfssSystemState(const size_t size):configuration_(size){
    configuration_.setToIdentity();
    configuration_.setPhaseToIdentity();
    parameters_.latticeSize=size;
  }

  explicit BfssSystemState(const BfssParameterSet& parameters):configuration_(parameters.latticeSize),parameters_(parameters){
    configuration_.setToIdentity();
    configuration_.setPhaseToIdentity();
    updateParameter();
  }

  BfssSystemState(BfssSystemState&)=default;
  BfssSystemState& operator=(BfssSystemState&)=default;

  BfssParameterSet& parameters(){
    return parameters_;
  }
  const BfssParameterSet& parameters()const{
    return parameters_;
  }

  StoredDataSet& data(){
    return storedDataSet_;
  }
  const StoredDataSet& data()const{
    return storedDataSet_;
  }

  BfssConfiguration& config(){
    return configuration_;
  }

  const BfssConfiguration& config()const{
    return configuration_;
  }

  auto& configNum(){
    return parameters_.confignumber;
  }

  /*
   * Action parameters need to be synchronized with DiracOperator parameters.
   */
  void updateParameter(){
    parameters().diracParameter.mu=parameters().mu;
    parameters().diracParameter.fermmass=parameters().fermmass;
    parameters().diracParameter.improved=parameters().improvedaction;
    parameters().diracParameter.gauged=!parameters().ungauged;
    parameters().diracParameter.fermionbc=parameters().fermionbc;
    parameters().diracParameter.setDefaults(parameters().latticeSpacing());
  }

  BfssDiracOperator getDiracOperator(){
    return BfssDiracOperator(&config(),&(parameters().diracParameter));
  }

  BfssDiracOperator getDiracOperatorDagger(){
    BfssDiracOperator ret(&config(),&(parameters().diracParameter));
    ret.setDag(true);
    return ret;
  }



};

struct BfssSystemStateFactory{
    using State=BfssSystemState;
    State* getState(Parser::ConfigParser& parser){
      BfssParameterSet parameter;
      readStateParameter(parser,parameter);
      return new State(parameter);
    }

};




#endif /* SRC_BFSS_BFSSSYSTEMSTATE_H_ */
