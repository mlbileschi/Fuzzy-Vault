#include "../inc/VaultMethod.h"
#include "../inc/RunTest.h"
#include <vector>
using namespace std;


int X_SHIFT;
int X_MASK;
int Y_MASK;
int THETA_MASK;
int X_VALUES;
int Y_VALUES;
int THETA_VALUES;

VaultBF::VaultBF(){
	INV = false;
	X_SHIFT = THETA_BITS + Y_BITS;
	X_MASK = ((int)pow(2,to_float(X_BITS)) -1) << X_SHIFT;
	Y_MASK = ((int)pow(2,to_float(Y_BITS)) -1) << THETA_BITS;
	THETA_MASK = pow(2,to_float(THETA_BITS)) -1;

	X_VALUES = (int)pow(2, to_float(X_BITS));
	Y_VALUES = (int)pow(2, to_float(Y_BITS));
	THETA_VALUES = (int)pow(2, to_float(THETA_BITS));

	//cout << "X_SHIFT: " << X_SHIFT << "   X_MASK: " << X_MASK << "   Y_MASK: " << Y_MASK  << "   THETA_MASK: " << THETA_MASK << endl;
}

void VaultBF::getDefaultDistro(vector<float>* distro, float* distroIncrement){
	*distroIncrement = 1;

	for(int i=0; i<3; i++){
		distro->push_back(1.0);
	}
	distro->push_back(0.0);
}


string VaultBF::getMethodType(){
	return "bf";
}


std::set<ZZ> VaultBF::minutiae2ZZ(const vector<minutia> &hTst){
	std::set<ZZ> ans;
	int n = hTst.size();
	for(int i=0; i<n; i++){ 

//This code was used to get close to uniform distributions for each variable and may be needed again		
/*
//cout << z;
static int poop[64] = {0};
//poop[(hTst[i].m_nX >> 3)&63]++;
//poop[toFeild((float)hTst[i].x, 250, 52, 6)]++;
poop[toFeild((float)hTst[i].y, 250, 77, 6)]++;
//poop[toFeild((float)hTst[i].theta, 180, 100, 4)]++;
//poop[toFeild(d2, 35, 20, 4)]++;
//poop[toFeild(sigma0, 40, 0, 1)]++;
//poop[toFeild(sigma1, 28, 0, 1)]++;
//poop[sigma2]++;
//poop[(sigma0/22)%8]++;
cout << "poop:";
for(int j=0; j<64; j++){
cout <<  " " << poop[j];
}
cout << endl;
*/
		ZZ z = to_ZZ(toFeild((float)hTst[i].x, 250, 52, X_BITS));
		z = z << Y_BITS;
		z = z + to_ZZ(toFeild((float)hTst[i].y, 250, 77, Y_BITS));
		z = z << THETA_BITS;
		z = z + to_ZZ(toFeild((float)hTst[i].theta, 180, 100, THETA_BITS)); 


		
	//	folded fold;
	//	fold.z1 = to_ZZ_p(z);
		ans.insert(z);
		
//	cout << "(" << (*hTst).marr[i].m_nX << ", " << (*hTst).marr[i].m_nY << ", " << (*hTst).marr[i].m_nTheta << ")-(";
//	cout << ((*hTst).marr[i].m_nX/8)%64 << ", " << ((*hTst).marr[i].m_nY/8)%64 << ", " << ((*hTst).marr[i].m_nTheta/25)%16 << ")   ";

	}
	return ans;
}



float VaultBF::distance(const ZZ_p& z1, const ZZ_p& z2){

// old code
/*
	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));

	int dx, dy, dtheta;

	dtheta = abs((zz1&15) - (zz2&15)); // <<2 to weight it closer to the distances?
	zz1 = zz1 >> 4;
	zz2 = zz2 >> 4;

	dy = abs((zz1&63) - (zz2&63));
	zz1 = zz1 >> 6;
	zz2 = zz2 >> 6;

	dx = abs((zz1&63) - (zz2&63));

	return sqrt(pow(to_float(dx),2) + pow(to_float(dy),2) + pow(to_float(dtheta),2));
*/


	//int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));


/*	return sqrt(
		pow(to_float(((zz1&64512) - (zz2&64512)) >> 10), 2) + 
		pow(to_float(((zz1&1008) - (zz2&1008)) >> 4), 2) + 
		pow(to_float((zz1&15) - (zz2&15)), 2)
		);
*/
	int dTheta = abs(((zz1&THETA_MASK) - (zz2&THETA_MASK)));
	if(dTheta > 180){dTheta = 360 - dTheta;}

	return root[
		(((zz1&X_MASK) - (zz2&X_MASK)) >> X_SHIFT) * (((zz1&X_MASK) - (zz2&X_MASK)) >> X_SHIFT)  + 
		(((zz1&Y_MASK) - (zz2&Y_MASK)) >> THETA_BITS) * (((zz1&Y_MASK) - (zz2&Y_MASK)) >> THETA_BITS)    + 
		dTheta * dTheta
		];




	//float dist = sqrt(pow(to_float(dx),2) + pow(to_float(dy),2) + pow(to_float(dtheta),2));

	//static int counter[25] = {0};
	//counter[(int)(dist/10)]++;
	
	//cout << "\r";
	//for(int i=0; i<12; i++)
	//	cout << counter[i] << " ";
	

	//return dist;

	
/*

	float x1, y1, t1, x2, y2, t2, dx, dy, dt;

	t1 = to_float(rep(z1)%to_ZZ(16));
	t2 = to_float(rep(z2)%to_ZZ(16));

	y1 = to_float((rep(z1)/to_ZZ(16))%to_ZZ(64));
	y2 = to_float((rep(z2)/to_ZZ(16))%to_ZZ(64));

	x1 = to_float(rep(z1)/to_ZZ(1024));
	x2 = to_float(rep(z2)/to_ZZ(1024));

	dx = abs(x2-x1);
	dy = abs(y2-y1);
	dt = abs(t2-t1);
	
	return sqrt(pow(dx,2) + pow(dy,2) + pow(dt,2));
*/

}




std::set<int> VaultBF::minutiae2int(const vector<minutia> &hTst){
	std::set<int> ans;
	int n = hTst.size();
	for(int i=0; i<n; i++){ 

//This code was used to get close to uniform distributions for each variable and may be needed again		
/*
//cout << z;
static int poop[64] = {0};
//poop[(hTst[i].m_nX >> 3)&63]++;
//poop[toFeild((float)hTst[i].m_nX, 190, 60, 6)]++;
//poop[toFeild((float)hTst[i].m_nY, 190, 77, 6)]++;
poop[toFeild((float)hTst[i].m_nTheta, 180, 100, 4)]++;
//poop[toFeild(d2, 35, 20, 4)]++;
//poop[toFeild(sigma0, 40, 0, 1)]++;
//poop[toFeild(sigma1, 28, 0, 1)]++;
//poop[sigma2]++;
//poop[(sigma0/22)%8]++;
cout << "poop:";
for(int j=0; j<16; j++){
cout <<  " " << poop[j];
}
cout << endl;
*/
/*
		int z = toFeild((float)hTst[i].x, 190, 60, X_BITS);
		z = z << Y_BITS;
		z = z + to_ZZ(toFeild((float)hTst[i].y, 190, 60, Y_BITS));
		z = z << THETA_BITS;
		z = z + to_ZZ(toFeild((float)hTst[i].theta, 180, 60, THETA_BITS)); 
*/

		int z = abs(hTst[i].x);
		z = z << Y_BITS;
		z = z + abs(hTst[i].y);
		z = z << THETA_BITS;
		z = z + abs(hTst[i].theta); 

		ans.insert(z);
		

	}

	return ans;
}






float VaultBF::distance(const int& z1, const int& z2){

	int dTheta = abs((z1&THETA_MASK) - (z2&THETA_MASK));
	if(dTheta > 180){dTheta = 360 - dTheta;}

	return root[
		((((z1&X_MASK) - (z2&X_MASK)) >> X_SHIFT) * (((z1&X_MASK) - (z2&X_MASK)) >> X_SHIFT))  + 
		((((z1&Y_MASK) - (z2&Y_MASK)) >> THETA_BITS) * (((z1&Y_MASK) - (z2&Y_MASK)) >> THETA_BITS))    + 
		(dTheta * dTheta)
		];
}






// bool VaultBF::compFromZero(const folded& fold1, const folded& fold2){return true;}
