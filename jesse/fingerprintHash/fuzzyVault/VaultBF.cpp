#include "VaultMethod.h"
#include <vector>
using namespace std;

VaultBF::VaultBF(){
	INV = false;
	distroIncrement = .3;
	
	distro.push_back(1.0);
	distro.push_back(1.0);
	distro.push_back(1.0);
	distro.push_back(1.0);
	distro.push_back(0.8);
	distro.push_back(0.5);
	distro.push_back(0.25);
	distro.push_back(0.2);
	distro.push_back(0.17);
	distro.push_back(0.1);
	for(int i=0; i<7; i++){
		distro.push_back(0.0);
	}
}


vector<folded> VaultBF::minutiae2ZZ(const vector<Minutiae> &hTst){
	vector<folded> ans;
	int n = hTst.size();
	for(int i=0; i<n; i++){ 

		ZZ z = to_ZZ((hTst[i].m_nX >> 3)&63);
		z = z << 6;
		z = z + to_ZZ((hTst[i].m_nY >> 3)&63);
		z = z << 4;
		z = z + to_ZZ((hTst[i].m_nTheta/25)&15);
		
		folded fold;
		fold.z1 = to_ZZ_p(z);
		(ans).push_back(fold);
		
//	cout << "(" << (*hTst).marr[i].m_nX << ", " << (*hTst).marr[i].m_nY << ", " << (*hTst).marr[i].m_nTheta << ")-(";
//	cout << ((*hTst).marr[i].m_nX/8)%64 << ", " << ((*hTst).marr[i].m_nY/8)%64 << ", " << ((*hTst).marr[i].m_nTheta/25)%16 << ")   ";

	}
	return ans;
}

/*
bool VaultBF::compFold(folded fold1, folded fold2){
	return ((rep(fold1.z1) < rep(fold2.z1)) > 0);
}

bool VaultBF::comp(ZZ_p z1, ZZ_p z2){
	return true;
}
*/

//void VaultMethod::sorter(vector<folded> &foldedVector){
//	sort(foldedVector.begin(), foldedVector.end(), compFold);
//}

float VaultBF::distance(const ZZ_p& z1, const ZZ_p& z2){

	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));

	int dx, dy, dtheta;

	dx = abs(zz1&63 - zz2&63);
	zz1 = zz1 >> 6;
	zz2 = zz2 >> 6;

	dy = abs(zz1&63 - zz2&63);
	zz1 = zz1 >> 6;
	zz2 = zz2 >> 6;

	dtheta = abs(zz1&15 - zz2&15);

	return sqrt(pow((float)(dx),2) + pow(to_float(dy),2) + pow(to_float(dtheta),2));

	
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