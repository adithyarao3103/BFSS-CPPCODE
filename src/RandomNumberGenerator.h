/*
 * RandomNumberGenerator.h
 *
 *  Created on: 11 Dec 2012
 *      Author: Georg Bergner
 *
 *  Just a simple wrapper for random number generators. Taken from a different project on simulations of Yang-Mills theories.
 */

#ifndef RANDOMNUMBERGENERATOR_H_
#define RANDOMNUMBERGENERATOR_H_
#include "basicdef.h"

#include <omp.h>
// Random number generator -- now also part of std.
#include <boost/random/ranlux.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

/**
 * Thread safe implementation of the random number generator.
 */
class RandomNumberGeneratorThreaded {
private:
	typedef boost::mt19937 RandomGen;
	typedef boost::uniform_01<Real> Distribution;
	typedef boost::variate_generator<RandomGen, Distribution > Random;
	Random* rnd_;
	Random** prnd_;
	size_t num_;
	const size_t separation_;

public:
	RandomNumberGeneratorThreaded(const RandomNumberGeneratorThreaded&)=delete;
	RandomNumberGeneratorThreaded& operator=(const RandomNumberGeneratorThreaded&)=delete;

	RandomNumberGeneratorThreaded():rnd_(NULL),prnd_(NULL),num_(0),separation_(100000) {
		num_=omp_get_max_threads();
		const unsigned long int default_seed=time(0);
		rnd_=new Random(RandomGen(default_seed),Distribution());
		prnd_ = new Random*[num_];
		for (size_t i(num_); i--;) {
			prnd_[i] = new Random(RandomGen(default_seed+i*12345678910),
					boost::uniform_01<Real>());
		}
	}

	virtual ~RandomNumberGeneratorThreaded() {
		delete rnd_;
		for(size_t i(num_);i--;){
			delete prnd_[i];
		}
		delete[] prnd_;
	}

	Real operator()(){
		return((*prnd_[omp_get_thread_num()])());
	}

	Random& getGenerator(){
		return (*rnd_);
	}

	size_t size()const{
		return (num_);
	}

	Random& getGenerator(size_t i){
		return (*prnd_[i]);
	}



};

class RandomNumberGenerator {
private:
  typedef boost::mt19937 RandomGen;
  typedef boost::uniform_01<Real> Distribution;
  typedef boost::variate_generator<RandomGen, Distribution > Random;
  Random rnd_;

public:
  RandomNumberGenerator(const RandomNumberGenerator&)=delete;
  RandomNumberGenerator& operator=(const RandomNumberGenerator&)=delete;

  RandomNumberGenerator():rnd_(RandomGen(time(0)),boost::uniform_01<Real>()) {
  }


  Real operator()(){
    return rnd_();
  }


};

class RandomNumberGeneratorReference {
private:
	RandomNumberGenerator* r_;
public:
	RandomNumberGeneratorReference(RandomNumberGenerator* r):r_(r){}
	RandomNumberGeneratorReference(const RandomNumberGeneratorReference& rr):r_(rr.r_){}

	Real operator()(){
		return ((*r_)());
	}

};

inline void drawGaussianRandom(Real* in,const size_t size, RandomNumberGenerator& rnd) {
	for (size_t i(0); i < size; i += 2) {
		const Real pref(sqrt(fabs(-2.0 * log(1.0 - rnd()))));
		const Real r(1.0 - rnd());
		in[i] = pref * sin(2 * pi * r);
		if(i+1<size)
		  in[i + 1] = pref * cos(2 * pi * r);
	}
}

inline void drawGaussianRandom(std::vector<Real>& in, RandomNumberGenerator& rnd) {
  if(in.size()==0)
    return;
  drawGaussianRandom(&in[0],in.size(),rnd);
}

inline void drawU1Random(Real* in,const size_t size, RandomNumberGenerator& rnd) {
	for (size_t i(0); i < size; i += 2) {
		// const Real pref(sqrt(fabs(-2.0 * log(1.0 - rnd()))));
		const Real r(1.0 - rnd());
		in[i] = sin(2 * pi * r);
		if(i+1<size)
		  in[i + 1] = cos(2 * pi * r);
	}
}

inline void drawU1Random(std::vector<Real>& in, RandomNumberGenerator& rnd) {
  if(in.size()==0)
    return;
  drawU1Random(&in[0],in.size(),rnd);
}







#endif /* RANDOMNUMBERGENERATOR_H_ */
