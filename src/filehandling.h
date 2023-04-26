/*
 * filehandling.h
 *
 *  Created on: 12.04.2023
 *      Author: Georg Bergner
 */

#ifndef SRC_BFSS_FILEHANDLING_H_
#define SRC_BFSS_FILEHANDLING_H_

#include "basicdef.h"
#include "bfsssystemstate.h"

#include <string>
#include <sstream>
#include <cstdlib> // system (Would be easier to use C++17 <filesystem>!)
#include <sys/stat.h>

// see also boost::filesystem (create_directory(), exists())

inline std::string commonDirectoryName(const BfssSystemState& state){
    std::ostringstream dirname;
    dirname << "_NC" << std::setfill('0') << std::setw(2) << numcolor
        << "_NI" << std::setfill('0') << std::setw(2)
        << BfssConfiguration::elementnum << "_NT"
        << std::setfill('0') << std::setw(3)
        << state.config().size() << "_L" << std::setfill('0')
        << std::setw(5) << state.parameters().lambda<<"_MU"<< std::setw(5) << state.parameters().mu;
    std::string dir("run" + dirname.str());
    return dir;
}


/*
 * This solution works most likely only on Unix OS.
 */
inline void initializeDirectoryIfNotPresent(const std::string& dirn){
  struct stat st;
  if(stat(dirn.c_str(),&st) != 0){
    std::string command("mkdir " + dirn);
    const int dummy=std::system(command.c_str()); // dummy to avoid warning.
  }
}


#endif /* SRC_BFSS_FILEHANDLING_H_ */
