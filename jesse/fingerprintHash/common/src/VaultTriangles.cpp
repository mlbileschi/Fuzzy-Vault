#include "../inc/VaultMethod.h"

VaultTriangles::VaultTriangles(){
	INV = true;
}


void VaultTriangles::getDefaultDistro(vector<float>* distro, float* distroIncrement){
	*distroIncrement = 1;

	for(int i=0; i<4; i++){
		distro->push_back(1.0);
	}
	distro->push_back(0);
}


std::set<ZZ> VaultTriangles::minutiae2ZZ(const vector<Minutiae>& temp){
	std::set<ZZ> ans;

	int n = temp.size();
//	cout << (*temp).nCnt << "  " << n << " . ";
	//n=15;
	for(int p=0; p<n; p++){ 
		


		double min[2];
		int index[2];

		for(int j=0; j<2; j++){
			min[j] = 1000.0;
			index[j] = 0;
		}

		for(int i=0; i<n; i++){
			if(i==p){continue;}
			double dx, dy, tee;

			ZZ delta_x = to_ZZ(abs(temp[i].m_nX - temp[p].m_nX));
			ZZ delta_y = to_ZZ(abs(temp[i].m_nY - temp[p].m_nY));
			//ZZ delta_theta = to_ZZ(abs(temp->marr[i].m_nTheta - temp->marr[p].m_nTheta));

			conv(dx, delta_x);
			conv(dy, delta_y);
			//conv(tee, delta_theta);

			double two_norm = sqrt(pow(dx,2) + pow(dy,2));

			if(two_norm > .01  &&  two_norm < min[1]){
				min[1] = two_norm;
				index[1] = i;

				for(int j=1; j>=0; j--){
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
		}


		float d0 = min[0];
		double d1 = min[1];


		double dx, dy, tee;

		ZZ delta_x = to_ZZ(abs(temp[index[0]].m_nX - temp[index[1]].m_nX));
		ZZ delta_y = to_ZZ(abs(temp[index[0]].m_nY - temp[index[1]].m_nY));

		conv(dx, delta_x);
		conv(dy, delta_y);

		double d2 = sqrt(pow(dx,2) + pow(dy,2));

		int sigma0 = abs(temp[index[0]].m_nTheta - temp[p].m_nTheta);
		int sigma1 = abs(temp[index[1]].m_nTheta - temp[p].m_nTheta);
		int sigma2 = abs(temp[index[0]].m_nTheta - temp[index[1]].m_nTheta);

		if(sigma0 > 180){sigma0 = 360 - sigma0;}
		if(sigma1 > 180){sigma1 = 360 - sigma1;}
		if(sigma2 > 180){sigma2 = 360 - sigma2;}
		sigma2 = abs(sigma0 - sigma1) == sigma2 ? 0 : 1;
			

//This code was used to get close to uniform distributions for each variable and may be needed again		
/*
//cout << z;
static int poop[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//poop[toFeild(d0, 25, 12, 4)]++;
//poop[toFeild(d1, 32, 14, 4)]++;
//poop[toFeild(d2, 35, 20, 4)]++;
//poop[toFeild(sigma0, 40, 0, 1)]++;
//poop[toFeild(sigma1, 28, 0, 1)]++;
//poop[sigma2]++;
//poop[(sigma0/22)%8]++;
cout << "poop:";
for(int i=0; i<16; i++){
cout <<  " " << poop[i];
}
cout << "   ";
*/
		ZZ z = to_ZZ(toFeild(d0, 25, 12, 4));
		z = z << 4;

		//z += to_ZZ((d1-32)/8)%to_ZZ(8);
		z += to_ZZ(toFeild(d1, 32, 14, 4));
		z = z << 4;

		//z = z + to_ZZ((d2)/10)%to_ZZ(8);
		z += to_ZZ(toFeild(d2, 35, 20, 4));
		z = z << 1;

		//z = z + to_ZZ((sigma0/22)%8);
		z += to_ZZ(toFeild(sigma0, 40, 0, 1));
		z = z << 1;

		//z = z + to_ZZ((sigma1)/22%8);
		z += to_ZZ(toFeild(sigma1, 28, 0, 1));
		z = z << 1;

		z += to_ZZ(sigma2);
		z = z << 1;

/*
		ZZ z = to_ZZ((d0-23)/8)%to_ZZ(8);
		z = z << 3;

		z += to_ZZ((d1-32)/8)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((d2)/10)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((sigma0/22)%8);
		z = z << 3;

		z = z + to_ZZ((sigma1)/22%8);
		z = z << 1;

		z += to_ZZ(sigma2);
*/

		//cout << z << " was z ";

		//z = z + to_ZZ(((*hTst).marr[i].m_nY/8)%64);
		//z = z << 6;
		//z = z + to_ZZ(((*hTst).marr[i].m_nTheta/25)%16);

	//	folded f;
	//	f.z1 = to_ZZ_p(z);
		ans.insert(z);
	}

	return ans;
}



float VaultTriangles::distance(const ZZ_p& z1, const ZZ_p& z2){

	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));
	int dd0, dd1, dd2, dsigma0, dsigma1, dsigma2;

	zz1 >>= 1;
	zz2 >>= 1;


	dsigma2 = abs((zz1&1) - (zz2&1)) << 3;
	zz1 >>= 1;
	zz2 >>= 1;

	dsigma1 = abs((zz1&1) - (zz2&1)) << 3;
	zz1 >>= 1;
	zz2 >>= 1;

	dsigma0 = abs((zz1&1) - (zz2&1)) << 3;
	zz1 >>= 1;
	zz2 >>= 1;

/*
	zz1 >>= 1;
	zz2 >>= 1;

if(zz1&7 != zz2&7){
	return 1000;
}
*/
//	zz1 >>= 3;
//	zz2 >>= 3;

	dd2 = abs((zz1&15) - (zz2&15));
	zz1 >>= 4;
	zz2 >>= 4;

	dd1 = abs((zz1&15) - (zz2&15));
	zz1 >>= 4;
	zz2 >>= 4;

	dd0 = abs((zz1&15) - (zz2&15));

	return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) + pow(to_float(dsigma0),2) + pow(to_float(dsigma1),2) + pow(to_float(dsigma2),2));
	//return to_float(dd0) + to_float(dd1) + to_float(dd2) + to_float(dsigma0) + to_float(dsigma1) + to_float(dsigma2);
	//return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) );
}



