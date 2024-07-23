/*
 *
 *  Created on: 28.02.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_BFSSCONFIGWRITER_H_
#define SRC_BFSS_BFSSCONFIGWRITER_H_

#include "bfsssystemstate.h"
#include "filehandling.h"

#include <fstream>
#include <vector>



class ConfigWriter {
private:
  std::string filename_;
  std::string directoryname_;
public:
  ConfigWriter() :
      filename_() {
  }
  void initialize(const BfssSystemState& state) {
    directoryname_=commonDirectoryName(state);
    filename_ = directoryname_ + "/config" + directoryname_;
    initializeDirectoryIfNotPresent(directoryname_);
  }


  void readConfig(BfssSystemState& state){
    readConfigurationFromOctaveFormat(state);
  }

  void writeConfig(const BfssSystemState& state){
    writeConfigurationToOctaveFormat(state);
  }

  void writeConfigurationToOctaveFormat(
      const BfssSystemState& state) {
    std::ostringstream cname;
    cname << "_" << std::setfill('0') << std::setw(8) << state.parameters().confignumber
        << "oct.txt";
    std::string thefilename(filename_ + cname.str());
    writeConfigurationToOctaveFormat(thefilename, state.config());
  }


  void writeConfigurationToOctaveFormat(
      const BfssSystemState& state,
      const std::vector<Real>& measurements,
      const std::list<std::string>& names) {
    std::ostringstream cname;
    cname << "_" << std::setfill('0') << std::setw(8) << state.parameters().confignumber
        << "oct.txt";
    std::string thefilename(filename_ + cname.str());
    writeConfigurationToOctaveFormat(thefilename, state.config());
    size_t ind(0);
    for (std::list<std::string>::const_iterator it(names.begin());
            it != names.end(); it++) {
      addMeasurementOctave(thefilename,measurements[ind++],*it);
    }

  }

  void readConfigurationFromOctaveFormat(BfssSystemState& state) {
    std::ostringstream cname;
    cname << "_" << std::setfill('0') << std::setw(8) << state.parameters().confignumber
        << "oct.txt";
    std::string thefilename((filename_ + cname.str()).c_str());
    readConfigurationFromOctaveFormat((filename_ + cname.str()),state.config());

  }



  void writeConfigurationToOctaveFormat(const std::string& filename,
      const BfssConfiguration& config) {
    std::ofstream file(filename.c_str());
    if (!file) {
      std::cout << "ERROR COULD NOT OPEN " << filename << std::endl;
      return;
    }
    file<<std::scientific<<std::setprecision(12);
    file << "# name: Phases\n" << "# type: matrix\n"
        << "# rows: 1\n" << "# columns: " << numcolor << std::endl;
        for (size_t jj = 0; jj < numcolor; ++jj) {
          file << config.phases()(jj)<<" ";
        }
    file << "\n";
    size_t size(config.size() * numcolor * numcolor);
    file << "# name: XField\n" << "# type: complex matrix\n" << "# rows: "
        << BfssConfiguration::elementnum << "\n"
        << "# columns: " << size << std::endl;
    for (size_t e = 0; e < config.elementnum; ++e) {
      for (size_t tt = 0; tt < config.size(); ++tt)
        for (size_t ii = 0; ii < numcolor; ++ii)
          for (size_t jj = 0; jj < numcolor; ++jj) {
            file << config(tt, e)(ii, jj);
          }
      file << "\n";
    }
    file << "\n";
  }


  void addMeasurementOctave(const std::string& filename,const Real value, const std::string& name){
    std::ofstream file(filename.c_str(),std::ios::out | std::ios::app);
    writeScalarOctave(file,name,value);
  }

  void addMeasurementOctave(const std::string& filename,const Complex value, const std::string& name){
    std::ofstream file(filename.c_str(),std::ios::out | std::ios::app);
    writeScalarOctave(file,name,value);
  }

    void writeScalarOctave(std::ostream& os,const std::string& name, const Real val) {
      os << "# name: " << name << "\n" << "# type: scalar\n" << " " << val << std::endl;
    }

    void writeScalarOctave (std::ostream& os, const std::string& name,const Complex val )
    {
      os<<"# name: "<<name<<"\n"<<"# type: complex scalar\n"<<" ("<<val.real() <<","<<val.imag() <<")"<<std::endl;
    }

  void readConfigurationFromOctaveFormat(const std::string& filename,
      BfssConfiguration& config) {
    std::ifstream file(filename.c_str());
    if (!file) {
      std::cout << "ERROR COULD NOT OPEN " << filename << std::endl;
      return;
    }

    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); //"# name: LinkConfig\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); //"# type: complex matrix\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# rows: 1\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# columns: "

    for (size_t jj = 0; jj < numcolor; ++jj) {
      file >> config.phases()(jj);
    }

    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# name: XField\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# type: complex matrix\n"
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# rows: "
    file.ignore(std::numeric_limits<std::streamsize>::max(),
        file.widen('\n')); // "# columns: "

    for (size_t e = 0; e < config.elementnum; ++e) {
      for (size_t tt = 0; tt < config.size(); ++tt)
        for (size_t ii = 0; ii < numcolor; ++ii)
          for (size_t jj = 0; jj < numcolor; ++jj) {
            file.ignore(10, file.widen('('));
            Real tmp;
            file >>tmp ;
            config(tt, e)(ii, jj).real(tmp);
            file.ignore(10, file.widen(','));
            file >> tmp;
            config(tt, e)(ii, jj).imag(tmp);
            file.ignore(10, file.widen(')'));
          }
    }
  }


};




#endif /* SRC_BFSS_BFSSCONFIGWRITER_H_ */
