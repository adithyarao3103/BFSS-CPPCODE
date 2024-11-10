// Code written by me to simulate SU(2) pure gauge theory in 1+1 dimensions. 
// This code is inefficient and prone to memory leaks. It is only to display an understanding of the lattice simulations
// Especially the definition of a lattice. Here I have implemented a more natural 4-dimensional array that mimics a physical lattice, while the most efficient way would be to encode the lattice serially.
// I have not implemented any observables so far. 
// The only part that has been implemented is the lattice action and the updating of the lattice according to the Markov Chain method.
// Have not used any standard libs in C++ for complex numbers and matrices. Implemented the necessary functionalities from scratch

#include <iostream>
#include <cmath>
#include <random>


#define nS 8
#define nT 8
#define threshold 0.5
#define N 2

//  Defining the complex number class

// Better implementation for all the classes would be to keep the class members private and only define the public functions. But it would take too much work to implement and I am short on time!

class complex {
public:
	double r, i;
	complex() {
		r = 0.0;
		i = 0.0;
	}
	complex(double real, double imag) {
		r = real;
		i = imag;
	}
	complex conjugate() {
		complex out;
		out.r = r;
		out.i = -1 * i;
		return (out);
	}
	void print() {
		std::cout << r << " + i" << i;
	}
};

complex operator*(complex c1, complex c2) {
		complex c;
		c.r = c1.r * c2.r - c1.i * c2.i;
		c.i = c1.r * c2.i + c1.i * c2.r;
		return(c);
	}

complex operator+(complex c1, complex c2) {
		complex c;
		c.r = c2.r + c1.r;
		c.i = c2.i + c1.i;
		return(c);
	}

complex operator*(double f, complex c) {
		complex out;
		out.r = c.r * f;
		out.i = c.i * f;
		return(out);
	}


// Defining the SU(2) matrix class

class su2 {
public:
	complex u[2][2];
	su2() {
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 2; col++) {
				u[row][col] = complex(0, 0);
			}
		}
	}
	su2(complex uin[2][2]) {
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 2; col++) {
				u[row][col] = uin[row][col];
			}
		}
	}
	su2 hermitian() {
		su2 out(*this);
		out.u[0][1] = u[1][0].conjugate();
		out.u[1][0] = u[0][1].conjugate();
		return(out);
	}
	complex trace() {
		return (u[0][0]+u[1][1]);
	}
	void print() {
		u[0][0].print();
		std::cout << "\t";
		u[0][1].print();
		std::cout << std::endl;
		u[1][0].print();
		std::cout << "\t";
		u[1][1].print();
		std::cout << std::endl;
	}
};

su2 operator+(su2 u1, su2 u2) {
		su2 sum;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				sum.u[i][j] = u1.u[i][j] + u2.u[i][j];
			}
		}
		return (sum);
	}

su2 operator*(complex c, su2 u) {
		su2 out;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				out.u[i][j] = c*u.u[i][j];
			}
		}
		return(out);
	}

su2 operator*(su2 u1, su2 u2) {
		su2 out;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				complex sum;
				for (int k = 0; k < 2; k++) {
					sum = sum + u1.u[i][k]*u2.u[k][j];
				}
				out.u[i][j] = sum;
			}
		}
		return(out);
	}

su2 operator*(double f, su2 u1) {
		su2 out;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				out.u[i][j] = f*u1.u[i][j];
			}
		}
		return(out);
	}

//  Some initializations

complex c_zero((double)0.0, (double)0.0);
complex c_one((double)1.0, (double)0.0);
complex c_minus_one((double)-1.0, (double)0.0);
complex c_i((double)0.0, (double)1.0);
complex c_minus_i = c_i.conjugate();

complex identity_[2][2] = {{c_one, c_zero}, {c_zero, c_one}};
complex sigx_[2][2] = {{c_zero, c_one}, {c_one, c_zero}};
complex sigy_[2][2] = {{c_zero, c_minus_i}, {c_i, c_zero}};
complex sigz_[2][2] = {{c_one, c_zero}, {c_zero, c_minus_one}};

su2 identity(identity_);
su2 sigx(sigx_);
su2 sigy(sigy_);
su2 sigz(sigz_);

std::random_device rd;
std::mt19937 gen(rd());
double LO = -0.5, HI = 0.5;
std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> distnS(0, nS-1);
std::uniform_int_distribution<std::mt19937::result_type> distnT(0, nT - 1);
std::uniform_int_distribution<std::mt19937::result_type> distmu(0, 1);

//  Defining the lattice class

class lattice {
public:
	su2*** lat = new su2**[nT];
	double beta;
	lattice(double beta) {
		this->beta = beta;

		//  Links unique to each site
		//	(t, s+1)
		//	   ^
		//	   |
		//	   | 2      
		//	   |      
		//	   |		0
		//	(t, s) ----------> (t+1, s)
		//  We also use the convention that the directional link connecting (t+1, s) to (t, s) is the hermitial conjugate of the directional link connecting (t, s) to (t+1, s).


		for (int t = 0; t < nT; t++) {
			lat[t] = new su2*[nS];
			for (int s = 0; s < nS; s++) {
				lat[t][s] = new su2[2];
				for (int mu = 0; mu < 2; mu++) {
					lat[t][s][mu] = identity;
				}
			}
		}		
	}
	void print() {
		for (int t = 0; t < nT; t++) {
			for (int s = 0; s < nS; s++) {
				for (int mu = 0; mu < 2; mu++) {
					std::cout << t << " " << s << " " << mu << std::endl;
					lat[t][s][mu].print();
				}
			}
		}
	}
	double action(){
		lattice l(*this);
		double S = 0;
		complex intermediate;
		for (int t = 0; t < nT; t++) {
			for (int s = 0; s < nS; s++) {
				su2 umunu = l.lat[t][s][0];
				// Below we use periodic boundary conditions in both spatial and temporal directions.
				if ( t == nT - 1 ) umunu = umunu*l.lat[0][s][1]; else umunu = umunu*l.lat[t + 1][s][1];
				if ( s == nS - 1 ) umunu = umunu*l.lat[t][0][0].hermitian(); else umunu = umunu*l.lat[t][s + 1][0].hermitian();
				umunu = umunu*l.lat[t][s][1].hermitian();
				intermediate = (identity + (-1)*umunu).trace();
				S += intermediate.r;
			}
		}
		S = S * beta / N;
		return (S);
	}
	void MCUpdate(){
		lattice lat1(*this);
		lattice lat2(*this);
		double initAction = lat1.action();
		double r0 = (double)LO + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (HI - LO)));
		double r1 = (double)LO + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (HI - LO)));
		double r2 = (double)LO + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (HI - LO)));
		double r3 = (double)LO + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (HI - LO)));
		double r = pow(pow(r1, 2) + pow(r2, 2) + pow(r3, 2), 0.5);
		double x0 = (r0 > 0 ? 1 : -1) * 0.8660254038;
		double x1 = 0.5 * r1 / r;
		double x2 = 0.5 * r2 / r;
		double x3 = 0.5 * r3 / r;
		su2 X = (x0*identity)+(x1*sigx)+(x2*sigy)+(x3*sigz);  // Gives a random SU(2) matrix which we will use to generate a new lattice from the given one 
		int temp = distnT(rng);
		int spat = distnS(rng);
		int direction = distmu(rng);
		lat2.lat[temp][spat][direction] = X*lat2.lat[temp][spat][direction];
		double finAction = lat2.action();
		double rndm = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double expDeltaS = exp(-(finAction - initAction));
		if ( rndm <= expDeltaS) // This is the markov chain selection criteria. If r <= exp(-\Delta S), where r \in [0,1) is a random number, the new lattice is accepted. Else, the old lattice is returned.
		{
			*this = lat2;
		}
	}
};

double Observable(lattice lat) {
	//  Here we implement the observable on the lattice.
	return (double)0.0;
}


int main()
{
	double beta=1;
	double actions[20];
	lattice lat(beta);
	int i = 0;
	while (i < 10001) {
		lat.MCUpdate();
		if (i % 100 == 0) std::cout <<"Count: "<<i <<" Action: "<<lat.action()<<" Observable: "<< Observable(lat) <<std::endl;
		i++;
	}
	//  The above code simply thermalises the lattice and prints the action. Any observable to be implemented must be done as a function of the lattice and should be run in a loop here. 
}



