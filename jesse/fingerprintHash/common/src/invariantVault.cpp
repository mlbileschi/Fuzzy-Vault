
#include "vault.h"
#include "invariantVault.h"
//#include "../../FingerprintHashing/random_gaussian.h"


bool comparer(triangles f1, triangles f2){return (f1.d[0] < f2.d[0]);}
bool compdouble(double f1, double f2){return f1<f2;}
double round(double d){return floor(d + 0.5);}
double truncate(double d){return round(d/TRUNK)*TRUNK;}

triangles Minutiae2triangles(int p, ptrFPTemplate &temp){

	triangles trist;

	double min[K];
	int index[K];

	for(int j=0; j<K; j++){
		min[j] = 1000.0;
		index[j] = 0;
	}

	for(int i=0; i<(*temp).nCnt; i++){
		if(i==p){continue;}
		double dx, dy, tee;

		ZZ delta_x = to_ZZ(abs(temp->marr[i].m_nX - temp->marr[p].m_nX));
		ZZ delta_y = to_ZZ(abs(temp->marr[i].m_nY - temp->marr[p].m_nY));
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
	}

	for(int j=0; j<K; j++){
		trist.d[j] = truncate(min[j]);
		trist.sigma[j] = (abs(temp->marr[index[j]].m_nTheta - temp->marr[p].m_nTheta));
		if(trist.sigma[j] > 180){trist.sigma[j] = 360 - trist.sigma[j];}
	}
  
	return trist;
}


ZZ_p Triangles2zzp(triangles tri){
	ZZ D = to_ZZ(0);
	ZZ S = to_ZZ(0);

	for(int u=0; u<K; u++){
		D += to_ZZ(ceil(tri.d[u]));
		S += to_ZZ(ceil(tri.sigma[u]));
	}

	ZZ z = D;
	z = z << 7;
	z = z + S;
	z = z << 1;
	return to_ZZ_p(z);
}


//Evaluates secret as a polynimial at z
ZZ_p function(ZZ_p z, vec_ZZ_p secret){
	ZZ_p sol=to_ZZ_p(0);
	ZZ_p zpow = to_ZZ_p(1);
	for(int i=0; i<secret.length(); i++){
		sol = sol + to_ZZ_p(rep(secret[i])*rep(zpow));
		zpow = zpow * z;
	}
	return sol;
}




void lockInv(vector<triangles>* vault, ptrFPTemplate hRef, vec_ZZ_p secret){
//cout << "\n________________________________________\n";
	//Add genuine points
	int n = hRef->nCnt;
	for(int i=0; i<n; i++){
		triangles trist;
		trist = Minutiae2triangles(i, hRef);
	
    	trist.chaff = false;

		ZZ_p z = Triangles2zzp(trist);

		trist.f1 = function(z, secret);
		if(FOLDED){trist.f2 = function(z + to_ZZ_p(1), secret);}


		vault->push_back(trist);

		//printf("X: %i\n",hRef->marr[i].m_nX);
		//printf("Y: %i\n",hRef->marr[i].m_nY);
		//printf("Theta: %i\n",hRef->marr[i].m_nTheta);

	}

	//Add chaff points
	double sum[K];
	double meandiff[K];
	double sig[K];
	for(int u=0; u<K; u++){
		sum[u] = 0;
		sig[u] = 0;
	}
	for(int y=0; y<vault->size(); y++){
		for(int r=0; r<K; r++){
			sum[r] += (*vault)[y].d[r];
			sig[r] += (*vault)[y].sigma[r];
		}
	}
	for(int u=0; u<K; u++){
		sum[u] = sum[u]/(double)vault->size();
		sig[u] = sig[u]/(double)vault->size();
	}


	meandiff[0] = sum[0];
	for(int u=1; u<K; u++){
		meandiff[u] = sum[u] - sum[u-1];
	}
	double devs[] = {13.93259, 16.91185, 21.47738, 22.89506, 22.89686, 26.23002, 28.57482, 28.98788, 30.04722, 30.37887, 31.56615, 31.91444, 31.74082, 32.20194, 31.6849, 32.03213, 32.02721, 31.68981, 31.87907, 32.30248};
	for(int i=0; i<NUM_CHAFF_POINTS; i++){
		triangles chaff;
		double disto = 0;
		double mean = 0;
		double sigo = 0;
		for(int y=0; y<K; y++){
			mean += meandiff[y];
		//	chaff.d[y] = truncate(randn_notrig(mean, devs[y]));
			chaff.d[y] = 10;

			if(rand()%3){
			//	sigo = (double)floor(randn_notrig(20, 10));
				sigo = 10;
			}else{
			//	sigo = (double)ceil(randn_notrig(160, 10));
				sigo = 10;
			}
			chaff.sigma[y] = sigo;

			//disto += (((double)rand())/((double)RAND_MAX))*(diff[y]*2);
			//chaff.d[y] = disto;
			//chaff.sigma[y] = (double) floor( (((double)rand())/((double)RAND_MAX))*180 );
			//cout << chaff.sigma[y] << " ";
		}
		//chaff.theta = to_ZZ_p(ceil((((double)rand())/((double)RAND_MAX))*360));
		sort(chaff.d, chaff.d+K);

		chaff.chaff = true;

		chaff.f1 = random_ZZ_p();

		ZZ_p zc = Triangles2zzp(chaff);
		if(chaff.f1 == function(zc, secret)){
			chaff.f1 = chaff.f1 + 300;
		}

		if(FOLDED){
			chaff.f2 = random_ZZ_p();
			if(chaff.f2 == function(zc + to_ZZ_p(1), secret) ){
				chaff.f2 = chaff.f2 + 300;
			}
		}


 vault->push_back(chaff);
	}

	//Sort using comp so genuine points aren't together
//****************	sort(vault->begin(), vault->end(), comparer);

}




/*******************************************************/

thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture unlockInv(vector<triangles>* vault, ptrFPTemplate hTst){

	thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture solution;
	solution.score = 0.0;
	int matches=0;
	//int score=0;
	vector<ZZ_p> x, fofx;
	vec_ZZ_p xvec, fofxvec;
	
	triangles test;

	double dee, es;

	int n = hTst->nCnt;
	int m = vault->size();



	int a, b, c, min_index, max_index;
	double maxx = 0.0;
	for(int i=0; i<n; i++){

		int tmp = 0;
		int maxkcount = 0;
		double minfar = 65537;
		int I = 0;

		test = Minutiae2triangles(i, hTst);

		for(int j=0; j<m; j++){

			int kcount = 0;
			double far = 0;
			bool fer[K];
			bool ger[K];
			for(int k=0; k<K; k++){
				fer[k] = false;
				ger[k] = false;
			}
			
			for(int k=0; k<K; k++){
				double mini = 65000.0;
				int indexi = 0;
				for(int t=0; t<K; t++){

					double delta_d = abs(test.d[k] - (*vault)[j].d[t]);
					double delta_s = abs(test.sigma[k] - (*vault)[j].sigma[t]);
				
					double two_norm = sqrt(delta_d*delta_d + delta_s*delta_s);
					if(two_norm < mini){mini = two_norm; indexi = t;}
				
				}
				if(mini < TOLERANCE){
					fer[indexi] = true;
					ger[k] = true;
					far = far + mini;
				}
			}
		
		int cfer = 0;
		int cger = 0;
			for(int d=0; d<K; d++){
				if(fer[d]){cfer++;}
				if(ger[d]){cger++;}
			}
			double farfer = 0;
			if(cfer < cger){kcount = cfer; farfer = cger;}else{kcount = cger; farfer = cfer;}

			if(kcount > maxkcount){
				maxkcount = kcount;
				minfar = far;
				I = j;
			}
			else if(kcount == maxkcount && (far/(1.0) < minfar)){
			//else if(kcount == maxkcount){
				minfar = far;
				I = j;
			}
		}//end j
	

		if(maxkcount >= KTHRESH){
		
				ZZ_p z2 = Triangles2zzp((*vault)[I]);
				//ZZ z = rep(z2);

				//check for duplicates
				int found = 0;
				
				for(int u=0; u<x.size(); u++){
					if (rep(x[u]) == rep(z2)){ found = 1;}
				}

				if(found == 0){
					x.push_back(z2);
					fofx.push_back((*vault)[I].f1);
					if(FOLDED){
						x.push_back(z2+to_ZZ_p(1));	
						fofx.push_back((*vault)[I].f2);
					}
					matches++;
					
					/*
					if((*vault)[I].chaff){
						solution.score = solution.score - Global_Parameters::DISTRO[floor(minL2)];}
					else{solution.score = solution.score + Global_Parameters::DISTRO[floor(minL2)];}
					*/

					if((*vault)[I].chaff){solution.score--;}else{solution.score++;}

				}
			
		}
	}//end i


	int num_points;

	if(FOLDED){
		num_points = 2*matches;
	}else{
		num_points = matches;
	}

	xvec.SetLength(num_points);
	fofxvec.SetLength(num_points);

	for(int i=0; i<num_points; i++){
		xvec[i] = x[i];
		fofxvec[i] = fofx[i];
	}

	ZZ_pX polyzz;

	solution.poly.SetLength(POLYNOMIAL_TERMS);
	int unlockable = POLYNOMIAL_TERMS;
	if(FOLDED){unlockable = unlockable/2;}

cout << "" << matches << "(" << solution.score << ") ";


		bool t = WelchBerlekamp (polyzz, xvec , fofxvec ,  (num_points-POLYNOMIAL_TERMS)/2);

		for(int i=0; i<POLYNOMIAL_TERMS; i++){
			solution.poly[i] = coeff(polyzz,i);
		}

	return solution;
	
}

