#include "../inc/VaultMethod.h"
#include "../inc/common.h"
#include <vector>
using namespace std;

VaultLetsTrySomething::VaultLetsTrySomething(){
	INV = true;

}


void VaultLetsTrySomething::getDefaultDistro(vector<float>* distro, float* distroIncrement){
	*distroIncrement = .3;
	
	distro->push_back(1.0);
	for(int i=0; i<7; i++){
		distro->push_back(0.0);
	}
}

std::set<ZZ> VaultLetsTrySomething::minutiae2ZZ(const vector<Minutiae> &hTst){

	int K = 8;
	int k = 8;

	std::set<ZZ> ans;

	double min[K];
	int index[K];

	for(int j=0; j<K; j++){
		min[j] = 1000.0;
		index[j] = 0;
	}
		
	ZZ delta_x, delta_y, delta_theta;
	double dx, dy, tee;
	for(int p=0; p<hTst.size(); p++){
		for(int i=0; i<hTst.size(); i++){
			if(i==p){continue;}

			ZZ delta_x = to_ZZ(abs(hTst[i].m_nX - hTst[p].m_nX));
			ZZ delta_y = to_ZZ(abs(hTst[i].m_nY - hTst[p].m_nY));
			//ZZ delta_theta = to_ZZ(abs(temp->marr[i].m_nTheta - temp->marr[p].m_nTheta));

			conv(dx, delta_x);
			conv(dy, delta_y);
			//conv(tee, delta_theta);

			double two_norm = sqrt(pow(dx,2) + pow(dy,2));

			if(two_norm > .01  &&  two_norm < min[K-1]){
				min[K-1] = two_norm;
				index[K-1] = i;

				for(int j=K-1; j>=0; j--){
					if(min[j] < min[j-1]){
						double tmpnorm = min[j-1];
						int tmpindex = index[j-1];
						min[j-1] = min[j];
						index[j-1] = index[j];
						min[j] = tmpnorm;
						index[j] = tmpindex;
					}
					else{
						break;
					}
				}
			}
		}//end i
		
		vector<double> two_normers;

		for(int i=0; i<k; i++){
			delta_x = to_ZZ(abs(hTst[index[i]].m_nX - hTst[p].m_nX));
			delta_y = to_ZZ(abs(hTst[index[i]].m_nY - hTst[p].m_nY));
			delta_theta = to_ZZ(abs(hTst[index[i]].m_nTheta - hTst[p].m_nTheta));

			conv(dx, delta_x);
			conv(dy, delta_y);
			conv(tee, delta_theta);

			two_normers.push_back(sqrt(pow(dx,2) + pow(dy,2) + pow(tee,2)));
		}

/*** This code was used to get close to uniform distributions for each variable and may be needed again
//cout << z;
static int poop[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
poop[toFeild(two_normers[0], 155, 80, 2)]++;
poop[toFeild(two_normers[1], 160, 80, 2)+4]++;
poop[toFeild(two_normers[2], 165, 80, 2)+8]++;
poop[toFeild(two_normers[3], 165, 85, 2)+12]++;
//poop[toFeild(two_normers[4], 165, 85, 2)]++;
//poop[toFeild(two_normers[5], 165, 85, 2)+4]++;
//poop[toFeild(two_normers[6], 165, 85, 2)+8]++;
//poop[toFeild(two_normers[7], 165, 85, 2)+12]++;
cout << "\npoop:";
for(int i=0; i<16; i++){
cout <<  " " << poop[i];
	if(i%4 == 3){cout << " -";}
}
cout << "   ";
*/

		ZZ z = to_ZZ(toFeild(two_normers[0], 155, 80, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[1], 160, 80, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[2], 165, 80, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[3], 165, 85, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[4], 165, 85, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[5], 165, 85, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[6], 165, 85, 2));
		z <<= 2;
		z += to_ZZ(toFeild(two_normers[7], 165, 85, 2));
		
		ans.insert(z);

	}// end p
	return ans;
}



float VaultLetsTrySomething::distance(const ZZ_p& z1, const ZZ_p& z2){

	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));

	int norm[8];

	norm[0] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[1] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[2] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[3] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[4] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[5] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[6] = abs((zz1&3) - (zz2&3)); 
	zz1 >>= 2;
	zz2 >>= 2;
	norm[7] = abs((zz1&3) - (zz2&3)); 

	int maxer = 0;
	for(int i=0; i<8; i++){
		maxer = norm[i] > maxer ? norm[i] : maxer;
	}
//cout << maxer << " ";
	return float(maxer);



}


