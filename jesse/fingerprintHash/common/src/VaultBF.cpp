#include "../inc/VaultMethod.h"
#include <vector>
using namespace std;

VaultBF::VaultBF(){
	INV = false;
}

void VaultBF::getDefaultDistro(vector<float>* distro, float* distroIncrement){
	*distroIncrement = 1;

	for(int i=0; i<3; i++){
		distro->push_back(1.0);
	}
	distro->push_back(0.0);
}


std::set<ZZ> VaultBF::minutiae2ZZ(const vector<Minutiae> &hTst){
	std::set<ZZ> ans;
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
cout << "\n";
*/

		ZZ z = to_ZZ(toFeild((float)hTst[i].m_nX, 190, 60, 6));
		z = z << 6;
		z = z + to_ZZ(toFeild((float)hTst[i].m_nY, 190, 60, 6));
		z = z << 4;
		z = z + to_ZZ(toFeild((float)hTst[i].m_nTheta, 180, 60, 4)); 
		
	//	folded fold;
	//	fold.z1 = to_ZZ_p(z);
		ans.insert(z);
		
//	cout << "(" << (*hTst).marr[i].m_nX << ", " << (*hTst).marr[i].m_nY << ", " << (*hTst).marr[i].m_nTheta << ")-(";
//	cout << ((*hTst).marr[i].m_nX/8)%64 << ", " << ((*hTst).marr[i].m_nY/8)%64 << ", " << ((*hTst).marr[i].m_nTheta/25)%16 << ")   ";

	}
	return ans;
}



float VaultBF::distance(const ZZ_p& z1, const ZZ_p& z2){

	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));

	int dx, dy, dtheta;

	dtheta = abs((zz1&15) - (zz2&15)) << 2; // <<2 to weight it closer to the distances?
	zz1 = zz1 >> 4;
	zz2 = zz2 >> 4;

	dy = abs((zz1&63) - (zz2&63));
	zz1 = zz1 >> 6;
	zz2 = zz2 >> 6;

	dx = abs((zz1&63) - (zz2&63));


	return sqrt(pow(to_float(dx),2) + pow(to_float(dy),2) + pow(to_float(dtheta),2));

	
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


// bool VaultBF::compFromZero(const folded& fold1, const folded& fold2){return true;}
