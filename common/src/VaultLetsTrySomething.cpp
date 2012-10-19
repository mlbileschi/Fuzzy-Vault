#include "../inc/VaultMethod.h"
#include "../inc/Common.h"
#include "../inc/RunTest.h"
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


string VaultLetsTrySomething::getMethodType(){
	return "new";
}


std::set<ZZ> VaultLetsTrySomething::minutiae2ZZ(const vector<minutia> &hTst){

	std::set<ZZ> ans;
	int n = hTst.size();
	float d0, d1, d2;
	
	for(int i=0; i<n; i++){
		for(int j=i+1; j<n; j++){
			for(int k=j+1; k<n; k++){

				d0 = root[((hTst[i].x - hTst[j].x) * (hTst[i].x - hTst[j].x)) + ((hTst[i].y - hTst[j].y) * (hTst[i].y - hTst[j].y))];
				d1 = root[((hTst[i].x - hTst[k].x) * (hTst[i].x - hTst[k].x)) + ((hTst[i].y - hTst[k].y) * (hTst[i].y - hTst[k].y))];
				d2 = root[((hTst[j].x - hTst[k].x) * (hTst[j].x - hTst[k].x)) + ((hTst[j].y - hTst[k].y) * (hTst[j].y - hTst[k].y))];

				if(d0 >100 || d1>100 || d2 >100){
					continue;
				}

/**** this code may be needed later
//cout << z;
static int poop[32] = {0};
//poop[toFeild(d0, 53, 23, 5)]++;
//poop[toFeild(d1, 53, 23, 5)]++;
poop[toFeild(d2, 53, 23, 5)]++;
cout << endl<< "poop:";
for(int i=0; i<32; i++){
cout <<  " " << poop[i];
	if(i%5 == 4){cout << " -";}
}
cout << "   ";
*/
	
				ZZ z = to_ZZ(toFeild(d0, 53, 23, 5));
				z = z << 5;

				//z += to_ZZ((d1-32)/8)%to_ZZ(8);
				z += to_ZZ(toFeild(d1, 53, 23, 5));
				z = z << 5;

				//z = z + to_ZZ((d2)/10)%to_ZZ(8);
				z += to_ZZ(toFeild(d2, 53, 23, 5));

				ans.insert(z);
}}} // end loops

/*
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

			ZZ delta_x = to_ZZ(abs(hTst[i].x - hTst[p].x));
			ZZ delta_y = to_ZZ(abs(hTst[i].y - hTst[p].y));
			//ZZ delta_theta = to_ZZ(abs(temp->marr[i].theta - temp->marr[p].theta));

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
			delta_x = to_ZZ(abs(hTst[index[i]].x - hTst[p].x));
			delta_y = to_ZZ(abs(hTst[index[i]].y - hTst[p].y));
			delta_theta = to_ZZ(abs(hTst[index[i]].theta - hTst[p].theta));

			conv(dx, delta_x);
			conv(dy, delta_y);
			conv(tee, delta_theta);

			two_normers.push_back(sqrt(pow(dx,2) + pow(dy,2) + pow(tee,2)));
		}
*/
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
cout << endl<< "poop:";
for(int i=0; i<16; i++){
cout <<  " " << poop[i];
	if(i%4 == 3){cout << " -";}
}
cout << "   ";
*/
/*
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
*/
	return ans;
}



float VaultLetsTrySomething::distance(const ZZ_p& z1, const ZZ_p& z2){


	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));
	int dd0, dd1, dd2, dsigma0, dsigma1, dsigma2;



/*
	zz1 >>= 1;
	zz2 >>= 1;

if(zz1&7 != zz2&7){
	return 1000;
}
*/
//	zz1 >>= 3;
//	zz2 >>= 3;

	dd2 = abs((zz1&31) - (zz2&31));
	zz1 >>= 5;
	zz2 >>= 5;

	dd1 = abs((zz1&31) - (zz2&31));
	zz1 >>= 5;
	zz2 >>= 5;

	dd0 = abs((zz1&31) - (zz2&31));

	return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2));


/*
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

*/

}



std::set<int> VaultLetsTrySomething::minutiae2int(const vector<minutia> &hTst){
	std::set<int> ans;
	return ans;
}

float VaultLetsTrySomething::distance(const int& z1, const int& z2){
	return 0;
}








