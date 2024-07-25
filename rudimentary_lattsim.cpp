// Code written by me for a simple simulation of SU(2) gauge theories on a 1+1 dimensional lattice. 
// There is a problem with the lattice not thermalizing but oscillating wildly and not reaching the ground state. 
// The action is reaching values like -9.94465e+012 and, on the other hand, 6.01714e+012. 
//I will take a deeper look at the problem when I have time!


#include <iostream>
#include <cmath>
#include <random>


#define nS 4
#define nT 4
#define threshold 0.5
#define beta 1
#define N 2

class complex {
public:
	float r, i;
	complex() {
		r = 0.0;
		i = 0.0;
	}
	complex(float real, float imag) {
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

complex operator*(float f, complex c) {
		complex out;
		out.r = c.r * f;
		out.i = c.i * f;
		return(out);
	}

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

su2 operator*(float f, su2 u1) {
		su2 out;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				out.u[i][j] = f*u1.u[i][j];
			}
		}
		return(out);
	}

complex c_zero((float)0.0, (float)0.0);
complex c_one((float)1.0, (float)0.0);
complex c_minus_one((float)-1.0, (float)0.0);
complex c_i((float)0.0, (float)1.0);
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
float LO = -0.5, HI = 0.5;
std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> distnS(0, nS-1);
std::uniform_int_distribution<std::mt19937::result_type> distnT(0, nT - 1);
std::uniform_int_distribution<std::mt19937::result_type> distmu(0, 3);


class lattice {
public:
	su2*** lat = new su2**[nT];
	lattice() {


		//		 1
		//	<----------
		//	|         ^
		//	| 3       | 2
		//	v         |
		//	---------->
		//	     0


		for (int t = 0; t < nT; t++) {
			lat[t] = new su2*[nS];

			for (int s = 0; s < nS; s++) {
				lat[t][s] = new su2[4];

				for (int mu = 0; mu < 4; mu++) {
					lat[t][s][mu] = identity;
				}
			}
		}		
	}

	void print() {
		for (int t = 0; t < nT; t++) {
			for (int s = 0; s < nS; s++) {
				for (int mu = 0; mu < 4; mu++) {
					std::cout << t << " " << s << " " << mu << std::endl;
					lat[t][s][mu].print();
				}
			}
		}
	}

};

float action(lattice l) {

		float S = 0;
		complex intermediate;


		for (int t = 1; t < nT - 1; t++) {

			for (int s = 1; s < nS - 1; s++) {

				su2 umunu;

				umunu = l.lat[t][s][0]*(l.lat[t + 1][s][2]*(l.lat[t + 1][s + 1][1]*l.lat[t][s + 1][3]));
				umunu = identity + umunu;

				intermediate = umunu.trace();

				S += intermediate.r;

			}

		}


		S = S * beta / N;

		return (S);

	}

lattice updateLattice(lattice lat2)
	{

		
		lattice origLat(lat2);

		float initAction = action(lat2);

		float r0 = (float)LO + static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX / (HI - LO)));
		float r1 = (float)LO + static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX / (HI - LO)));
		float r2 = (float)LO + static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX / (HI - LO)));
		float r3 = (float)LO + static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX / (HI - LO)));

		float r = pow(pow(r1, 2) + pow(r2, 2) + pow(r3, 2), 0.5);

		float x0 = (r0 > 0 ? 1 : -1) * 0.8660254038;
		float x1 = 0.5 * r1 / r;
		float x2 = 0.5 * r2 / r;
		float x3 = 0.5 * r3 / r;

		su2 X = (x0*identity)+(x1*sigx)+(x2*sigy)+(x3*sigz);

		int temp = distnT(rng);
		int spat = distnS(rng);
		int direction = distmu(rng);

		// std::cout<< temp << " " << spat << " " << direction << std::endl;
		// X.print();

		// lat2.lat[temp][spat][direction].print();

		lat2.lat[temp][spat][direction] = X*lat2.lat[temp][spat][direction];

		// lat2.lat[temp][spat][direction].print();

		if (direction == 0){
			if (temp == nT-1){
				lat2.lat[0][spat][1] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
			else{
				lat2.lat[temp+1][spat][1] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
		}
		else if (direction == 1){
			if (temp == 0){
				lat2.lat[nT-1][spat][0] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
			else{
				lat2.lat[temp-1][spat][0] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
		} else if (direction == 2){
			if (spat == nS-1){
				lat2.lat[temp][0][3] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
			else{
				lat2.lat[temp][spat+1][3] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
		}
		else if (direction == 3){
			if (spat == 0){
				lat2.lat[temp][nS-1][2] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
			else{
				lat2.lat[temp][spat-1][2] = (X*lat2.lat[temp][spat][direction]).hermitian();
			}
		}

		float finAction = action(lat2);

		float rndm = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float expDeltaS = exp(-(finAction - initAction));
		// std::cout<<"Random: "<<rndm<<" Final: "<<finAction<<" Initial: "<<initAction<<" Exp: "<<expDeltaS<<std::endl;

		if ( rndm <= expDeltaS)
		{
			// std::cout<<"Change"<<std::endl;
			return lat2;
		}
		else {
			// std::cout<<"No Change"<<std::endl;
			return origLat;
		}
	}

int main()
{
	// complex c1(10, 20);
	// complex c2(c1);
	// c2 = complex(1,2);
	// c2.print();
	// c1.print();
	lattice lat;
	int i = 0;
	while (i < 5000) {
		lat = updateLattice(lat);
		// std::cout <<i <<" "<< action(lat)<<std::endl;
		i++;
	}
		std::cout <<i <<" "<< action(lat)<<std::endl;

	// complex c1_[2][2] = {{c_one,2*c_one}, {3*c_one,4*c_one}};
	// complex c2_[2][2] = {{4*c_one, 5*c_one}, {6*c_one, 10*c_one}};

	// su2 c1(c1_), c2(c2_);

	// (c1*c2).print();
	
	
	// complex onetwo(1,2), twoone(2,1);
	// complex test_[2][2] = {{c_zero, onetwo}, {twoone, c_zero}};
	// su2 test(test_);
	// test.print();
	// test.hermitian();
	// test.print();

}





//class latticeNode {
//	su2 U_tp;
//	su2 U_xp;
//
//	latticeNode* xp;
//	latticeNode* tp;
//	latticeNode* xm;
//	latticeNode* tm;
//};
//
//class lattice {
//private:
//	latticeNode* origin = new latticeNode();
//public:
//	lattice() {
//		origin->U_tp = identity;
//		origin->U_xp = identity;
//
//		latticeNode* current = origin;
//		
//		/*
//		create a connected lattice like this 
//
//		x
//		o o o 
//		| | | 
//		o o o 
//		| | |
//		o-o-o t
//		
//		*/
//
//
//		for (t = 0; t < nT-1; t++) {
//
//			latticeNode* beginningX = current;
//
//			for (i = 0; i < nS-1; i++) {
//				current->xp = new latticeNode();
//				current->xp->xm = current;
//				current = current->xp;
//				current->U_xp = identity;
//				current->U_tp = identity;
//			}
//
//			current->xp = new latticeNode();
//			current = current->xp;
//			current->xp = beginningX;
//			beginningX->xm = current;
//			current->U_xp = identity;
//			current->U_tp = identity;
//			current = beginningX;
//
//
//			current->tp = new latticeNode();
//			current->tp->tm = current;
//			current = current->tp;
//			current->U_xp = identity;
//			current->U_tp = identity;
//
//		}
//		ccurrent->tp = new latticeNode();
//		current->tp->tm = current;
//		current = current->tp;
//		origin->tm = current;
//		current->U_xp = identity;
//		current->U_tp = identity;
//
//
//		//link the nodes in t direction
//
//		latticeNode* xNode = new latticeNode();
//		latticeNode* tNode = new latticeNode();
//		latticeNode* tNodePrev = new latticeNode();
//		latticeNode* tNodeNext = new latticeNode();
//		latticeNode* xNodePrev = new latticeNode();
//		latticeNode* xNodeNext = new latticeNode();
//
//		tNode = origin;
//		tNodePrev = tNode->tm;
//		tNodeNext = tNode->tp;
//
//
//		for (t = 0; t < nT; t++) {
//
//
//			xNode = tNode->xp;
//			xNodePrev = tNodePrev->xp;
//			xNodeNext = tNodeNext->xp;
//
//			for (i = 0; i < nS; i++) {
//			
//				xNode->tm = xNodePrev;
//				xNode->tp = xNodeNext;
//				xNodePrev = xNodePrev->xp;
//				xNode = xNode->xp;
//				xNodeNext = xNodeNext->xp;
//
//			}
//
//			tNodePrev = tNode;
//			tNode = tNodeNext;
//			tNodeNext = tNodeNext->tp;
//
//		}
//	}
//
//	void updateLattice() {
//		latticeNode* updateNode = origin;
//		for (i = 0; i < distnS(rng); i++) {
//			updateNode = updateNode->xp;
//		}
//		for (i = 0; i < distnT(rng); i++) {
//			updateNode = updateNode->tp;
//		}
//
//
//		float r0 = LO + static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
//		float r1 = LO + static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
//		float r2 = LO + static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
//		float r3 =  LO + static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
//
//		float r = pow(pow(r1, 2) + pow(r2, 2) + pow(r3, 2), 0.5);
//
//		float x0 = (r0 > 0 ? 1 : -1) * 0.8660254038;
//		float x1 = 0.5 * r1 / r;
//		float x2 = 0.5 * r2 / r;
//		float x3 = 0.5 * r3 / r;
//
//		su2 X = (identity.scalarMul(x0)).add((sigx.scalarMul(x1)).add((sigy.scalarMul(x2)).add((sigz.scalarMul(x3))));
//
//		if (dist4(rng) == 1) {
//			updateNode->tp = updateNode->tp.leftMul(X);
//		}
//
//	}
//
//	
//};


