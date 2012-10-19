#include "VaultMethod.h"


VaultTriangles::VaultTriangles(){
	INV = true;
	distroIncrement = 1.0;
	distro.push_back(1.0);
	//distro.push_back(0.7);
	//distro.push_back(0.4);
	//distro.push_back(0.15);
	for(int i=0; i<600; i++){
		distro.push_back(1.0);
	}
}

vector<folded> VaultTriangles::minutiae2ZZ(const vector<Minutiae>& temp){
	vector<folded> ans;

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


		double d0 = min[0];
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

		if(abs(sigma0 - sigma1) == sigma2){
			sigma2 = 0;
		}else{sigma2=1;}
		
//cout << 16.8%2 << " ";
	//	ZZ z = to_ZZ(d0);
		ZZ z = to_ZZ((d0-20)/10)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((d1-25)/10)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((d2)/20)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((sigma0)/45)%to_ZZ(8);
		z = z << 3;

		z = z + to_ZZ((sigma1)/45)%to_ZZ(8);
		z = z << 1;

		z = z + to_ZZ(sigma2);

		//cout << z << " was z ";

		//z = z + to_ZZ(((*hTst).marr[i].m_nY/8)%64);
		//z = z << 6;
		//z = z + to_ZZ(((*hTst).marr[i].m_nTheta/25)%16);

		folded f;
		f.z1 = to_ZZ_p(z);
		ans.push_back(f);
	}

	return ans;
}



float VaultTriangles::distance(const ZZ_p& z1, const ZZ_p& z2){

	int zz1, zz2;
	zz1 = to_int(rep(z1));
	zz2 = to_int(rep(z2));

	int dd0, dd1, dd2, dsigma0, dsigma1, dsigma2;

	dd0 = abs(zz1&7 - zz2&7);
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dd1 = abs(zz1&7 - zz2&7);
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dd2 = abs(zz1&7 - zz2&7);
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dsigma0 = abs(zz1&7 - zz2&7);
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dsigma1 = abs(zz1&7 - zz2&7);
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;
	
	dsigma2 = abs(zz1 - zz2);

/*
	dd1 = abs(zz1&56 - zz2&56);
	dd2 = abs(zz1&448 - zz2&448);
	dsigma0 = abs(zz1&3584 - zz2&3584);
	dsigma1 = abs(zz1&28672 - zz2&28672);
	dsigma2 = abs(zz1&32768 - zz2&32768);
	*/
//cout << dsigma2 << "";
	return sqrt(pow(to_float(dd0),2) + pow(to_float(dd1),2) + pow(to_float(dd2),2) + pow(to_float(dsigma0),2) + pow(to_float(dsigma1),2) + pow(to_float(dsigma2),2));
	//return sqrt(pow(dd0,2) + pow(dd1,2) + pow(dd2,2) + pow(dsigma0,2) + pow(dsigma1,2) + pow(dsigma2,2));
/*
	dd0 = abs( to_float(zz1%eight) - to_float(zz2%eight) );
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dd1 = abs( to_float(zz1%to_ZZ(8)) - to_float(zz2%to_ZZ(8)) );
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dd2 = abs( to_float(zz1%to_ZZ(8)) - to_float(zz2%to_ZZ(8)) );
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dsigma0 = abs( to_float(zz1%to_ZZ(8)) - to_float(zz2%to_ZZ(8)) );
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dsigma1 = abs( to_float(zz1%to_ZZ(8)) - to_float(zz2%to_ZZ(8)) );
	zz1 = zz1 >> 3;
	zz2 = zz2 >> 3;

	dsigma2 = abs( to_float(zz1%to_ZZ(2)) - to_float(zz2%to_ZZ(2)) );


	
	return sqrt(pow(dd0,2) + pow(dd1,2) + pow(dd2,2) + pow(dsigma0,2) + pow(dsigma1,2) + pow(dsigma2,2));
*/
}



