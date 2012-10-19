#include "../inc/VaultMethod.h"



int D0_BITS = 6;
int D1_BITS = 6;
int D2_BITS = 6;
int SIGMA0_BITS = 6;
int SIGMA1_BITS = 6;
int SIGMA2_BITS = 1;

int D0_SHIFT;
int D1_SHIFT;
int D2_SHIFT;
int SIGMA0_SHIFT;
int SIGMA1_SHIFT;
int SIGMA2_SHIFT;

int D0_MASK;
int D1_MASK;
int D2_MASK;
int SIGMA0_MASK;
int SIGMA1_MASK;
int SIGMA2_MASK;


VaultTriangles::VaultTriangles(){

	INV = true;
	D0_SHIFT = D1_BITS + D2_BITS + SIGMA0_BITS + SIGMA1_BITS + SIGMA2_BITS;
	D1_SHIFT = D2_BITS + SIGMA0_BITS + SIGMA1_BITS + SIGMA2_BITS;
	D2_SHIFT = SIGMA0_BITS + SIGMA1_BITS + SIGMA2_BITS;
	SIGMA0_SHIFT = SIGMA1_BITS + SIGMA2_BITS;
	SIGMA1_SHIFT = SIGMA2_BITS;
	SIGMA2_SHIFT = 0;

	D0_MASK = ((int)pow(2,to_float(D0_BITS)) -1) << D0_SHIFT;
	D1_MASK = ((int)pow(2,to_float(D1_BITS)) -1) << D1_SHIFT;
	D2_MASK = ((int)pow(2,to_float(D2_BITS)) -1) << D2_SHIFT;
	SIGMA0_MASK = ((int)pow(2,to_float(SIGMA0_BITS)) -1) << SIGMA0_SHIFT;
	SIGMA1_MASK = ((int)pow(2,to_float(SIGMA1_BITS)) -1) << SIGMA1_SHIFT;
	SIGMA2_MASK = ((int)pow(2,to_float(SIGMA2_BITS)) -1) << SIGMA2_SHIFT;

}


void VaultTriangles::getDefaultDistro(vector<float>* distro, float* distroIncrement){
	*distroIncrement = 1;

	for(int i=0; i<4; i++){
		distro->push_back(1.0);
	}
	distro->push_back(0);
}


string VaultTriangles::getMethodType(){
	return "triangles";
}


std::set<ZZ> VaultTriangles::minutiae2ZZ(const vector<minutia>& temp){
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

			ZZ delta_x = to_ZZ(abs((temp[i]).x - (temp[p]).x));
			ZZ delta_y = to_ZZ(abs(temp[i].y - temp[p].y));
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

		ZZ delta_x = to_ZZ(abs(temp[index[0]].x - temp[index[1]].x));
		ZZ delta_y = to_ZZ(abs(temp[index[0]].y - temp[index[1]].y));

		conv(dx, delta_x);
		conv(dy, delta_y);

		double d2 = sqrt(pow(dx,2) + pow(dy,2));

		int sigma0 = abs(temp[index[0]].theta - temp[p].theta);
		int sigma1 = abs(temp[index[1]].theta - temp[p].theta);
		int sigma2 = abs(temp[index[0]].theta - temp[index[1]].theta);

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
		ZZ z = to_ZZ(toFeild(d0, 25, 12, D0_BITS));
		z = z << D1_BITS;

		//z += to_ZZ((d1-32)/8)%to_ZZ(8);
		z += to_ZZ(toFeild(d1, 32, 14, D1_BITS));
		z = z << D2_BITS;

		//z = z + to_ZZ((d2)/10)%to_ZZ(8);
		z += to_ZZ(toFeild(d2, 35, 20, D2_BITS));
		z = z << SIGMA0_BITS;

		//z = z + to_ZZ((sigma0/22)%8);
		z += to_ZZ(toFeild(sigma0, 40, 0, SIGMA0_BITS));
		z = z << SIGMA1_BITS;

		//z = z + to_ZZ((sigma1)/22%8);
		z += to_ZZ(toFeild(sigma1, 28, 0, SIGMA0_BITS));
		z = z << SIGMA2_BITS;

		z += to_ZZ(sigma2);
		//z = z << 1;

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


	dsigma2 = abs((zz1&SIGMA2_MASK) - (zz2&SIGMA2_MASK)) >> SIGMA2_SHIFT;
	dsigma1 = abs((zz1&SIGMA1_MASK) - (zz2&SIGMA1_MASK)) >> SIGMA1_SHIFT;
	dsigma0 = abs((zz1&SIGMA0_MASK) - (zz2&SIGMA0_MASK)) >> SIGMA0_SHIFT;
	dd2 = abs((zz1&D2_MASK) - (zz2&D2_MASK)) >> D2_SHIFT;
	dd1 = abs((zz1&D1_MASK) - (zz2&D1_MASK)) >> D1_SHIFT;
	dd0 = abs((zz1&D0_MASK) - (zz2&D0_MASK)) >> D0_SHIFT;

	return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) + pow(to_float(dsigma0),2) + pow(to_float(dsigma1),2) + pow(to_float(dsigma2),2));

	//return to_float(dd0) + to_float(dd1) + to_float(dd2) + to_float(dsigma0) + to_float(dsigma1) + to_float(dsigma2);
	//return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) );
}






std::set<int> VaultTriangles::minutiae2int(const vector<minutia> &hTst){
	std::set<int> ans;

	int n = hTst.size();
//	cout << (*hTst).nCnt << "  " << n << " . ";
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

			int delta_x = abs((hTst[i]).x - (hTst[p]).x);
			int delta_y = abs(hTst[i].y - hTst[p].y);
			//int delta_theta = abs(hTst->marr[i].m_nTheta - hTst->marr[p].m_nTheta);

			//conv(dx, delta_x);
			//conv(dy, delta_y);
			//conv(tee, delta_theta);

			double two_norm = sqrt(pow((double)delta_x,2.0) + pow((double)delta_y,2.0));

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

		int delta_x = abs(hTst[index[0]].x - hTst[index[1]].x);
		int delta_y = abs(hTst[index[0]].y - hTst[index[1]].y);

		//conv(dx, delta_x);
		//conv(dy, delta_y);

		double d2 = sqrt(pow((double)delta_x,2.0) + pow((double)delta_y,2.0));

		int sigma0 = abs(hTst[index[0]].theta - hTst[p].theta);
		int sigma1 = abs(hTst[index[1]].theta - hTst[p].theta);
		int sigma2 = abs(hTst[index[0]].theta - hTst[index[1]].theta);

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
		int z = toFeild(d0, 25, 12, D0_BITS);
		z = z << D1_BITS;

		z += toFeild(d1, 32, 14, D1_BITS);
		z = z << D2_BITS;

		z += toFeild(d2, 35, 20, D2_BITS);
		z = z << SIGMA0_BITS;

		z += toFeild(sigma0, 40, 0, SIGMA0_BITS);
		z = z << SIGMA1_BITS;

		z += toFeild(sigma1, 28, 0, SIGMA0_BITS);
		z = z << SIGMA2_BITS;

		z += sigma2;


		ans.insert(z);
	}

	return ans;
}






float VaultTriangles::distance(const int& z1, const int& z2){

	int zz1, zz2;
	zz1 = z1;
	zz2 = z2;
	int dd0, dd1, dd2, dsigma0, dsigma1, dsigma2;


	dsigma2 = abs((zz1&SIGMA2_MASK) - (zz2&SIGMA2_MASK)) >> SIGMA2_SHIFT;
	dsigma1 = abs((zz1&SIGMA1_MASK) - (zz2&SIGMA1_MASK)) >> SIGMA1_SHIFT;
	dsigma0 = abs((zz1&SIGMA0_MASK) - (zz2&SIGMA0_MASK)) >> SIGMA0_SHIFT;
	dd2 = abs((zz1&D2_MASK) - (zz2&D2_MASK)) >> D2_SHIFT;
	dd1 = abs((zz1&D1_MASK) - (zz2&D1_MASK)) >> D1_SHIFT;
	dd0 = abs((zz1&D0_MASK) - (zz2&D0_MASK)) >> D0_SHIFT;

	return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) + pow(to_float(dsigma0),2) + pow(to_float(dsigma1),2) + pow(to_float(dsigma2),2));

}






