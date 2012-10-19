
#include "vault.h"
//#include <set>


ZZ_p Vault::function(ZZ_p z, vec_ZZ_p secret){
	ZZ_p sol=to_ZZ_p(0);
	ZZ_p zpow = to_ZZ_p(1);
	for(int i=0; i<secret.length(); i++){
		sol = sol + to_ZZ_p(rep(secret[i])*rep(zpow));
		zpow = zpow * z;
	}
	return sol;
}

//


Vault::Vault(VaultMethod* m){
	method = m;
}

void Vault::lock(vector<Minutiae> hRef, vec_ZZ_p secret){

	//int n = hRef.size();
	//Add genuine points
//	cout << (*hRef).nCnt << "  " << n << " . ";
	vector<folded> zValues = method->minutiae2ZZ(hRef);
	//VaultMethod method = VaultBF();
	folded fold;
	fold.chaff = false;		
//	for(int i=0; i<n; i++){
	for(int i=0; i<(zValues).size(); i++){
		fold.z1 = (zValues)[i].z1;
		fold.f1 = function((zValues)[i].z1, secret);

		vault.push_back(fold);
	}


	//Add chaff points

	for(int i=0; i<NUM_CHAFF_POINTS; i++){

		folded chaff;
		chaff.chaff = true;

		chaff.z1 = random_ZZ_p();
		chaff.f1 = random_ZZ_p();

		if(chaff.f1 == function(chaff.z1, secret)){
			chaff.f1 = chaff.f1 + 2; //arbitrary
		}
/*
		if(FOLDED){
			chaff.f2 = random_ZZ_p();
			if(chaff.f2 == function(to_ZZ_p(zc+to_ZZ(1)), secret)){
				chaff.f2 = chaff.f2 + 2;
			}
		}
*/

		vault.push_back(chaff);
	}

	//Sort using comp so genuine points aren't together

	Compare_functor free_compare(method);
		sort(vault.begin(), vault.end(), free_compare); 
/*
	for(int q=0; q<vault.size(); q++){
		cout << (vault)[q].z1 << "\n";
	}
	cout << "\n";
*/
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////

thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture Vault::unlock(vector<Minutiae> hTst){

	thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture data;
	data.score = 0.0;

	vector<folded> candidate = (method->minutiae2ZZ(hTst));

	Compare_functor free_compare(method);
	sort(candidate.begin(), candidate.end(), free_compare);


	//sort(candidate.begin(), candidate.end(), this()); 
	//cout << " !!!" << method->sorter(5) << "!!! ";
	//method->sorter(candidate);
//	for(int q=0; q<candidate.size(); q++){
//		//cout << (candidate)[q].z1 << "\n";
//		ZZ z1  = rep((candidate)[q].z1);
//		int ctr = 0;
//		for ( long long int mask = 1<<15; mask>=1; mask=mask/2){
//			if(ctr%6==0) cout<<" ";
//			((z1 & mask)!=0) ? cout<<1 : cout<<0;
//			ctr++;
//		}
//		cout<<endl;
//	}
////	cout << "\n";

	//compareExtract
	vector<float> distro = method->getDistro();

	int maxDistro = distro.size();
	float j = 0;
	while(j==0){
		maxDistro--;
		j = distro[maxDistro];
	}

	maxDistro = maxDistro*method->getDistroIncrement();

	std::set<int> indexPoints;

	int base = 0;

	int n = (vault).size();
	int m = candidate.size();

	float canN;
	float vaultN;
	float dist;

	float minDist;
	int minIndex;
	vector<float> vaultNs;

	for(int g=0; g<n; g++){
		vaultNs.push_back(method->distance((vault)[g].z1, to_ZZ_p(0)));
	}

	for(int i=0; i<m; i++){
		canN = method->distance((candidate)[i].z1, to_ZZ_p(0));
		minDist = maxDistro+method->getDistroIncrement(); // make it a little bigger than the max
		minIndex = -1;
		for(int k=base; k<n; k++){
			//vaultN = method->distance((vault)[k].z1, to_ZZ_p(0)); 
			if(canN - maxDistro > vaultNs[k]){
				base++;
				continue;
			}
			if(canN + maxDistro < vaultNs[k]){
				break;
			}
			dist = method->distance((candidate)[i].z1, (vault)[k].z1);
			if(dist<minDist){
				minDist = dist;
				//cout << minDist << " ";
				minIndex = k;
			}
		}//end k
//		cout << minDist << " ";
		//TODO: random by distro 
		if(minIndex != -1){
			float rando = (float)rand()/(float)RAND_MAX;
			if(rando < distro[(int)ceil(minDist/method->getDistroIncrement())]){
				//cout << "    " << rando << " " << minDist << " " << distro[(int)ceil(minDist/method->getDistroIncrement())];
				//cout << " " << (int)ceil(minDist/method->getDistroIncrement()) << "   ";
				indexPoints.insert(minIndex);
			}
		}
	}//end i

	//@TODO score by distro

	int numPoints = indexPoints.size();

	vec_ZZ_p xvec, fofxvec;
	xvec.SetLength(numPoints);
	fofxvec.SetLength(numPoints);

	int numgen=0, numimp=0;
	std::set<int>::iterator it;
	int i=0;
	for(it = indexPoints.begin(); it != indexPoints.end(); it++ ) {
		xvec[i] = (vault)[*it].z1;
		fofxvec[i] = (vault)[*it].f1;
		i++;
		if((vault)[*it].chaff){
			numimp++;
			data.score--;
		}else{
			numgen++;
			data.score++;
		}
	}  

	ZZ_pX polyzz;

	data.poly.SetLength(POLYNOMIAL_TERMS);
	if((numPoints-POLYNOMIAL_TERMS)/2 > 0){
		bool t = WelchBerlekamp (polyzz, xvec , fofxvec ,  (numPoints-POLYNOMIAL_TERMS)/2);
	}

		for(int i=0; i<POLYNOMIAL_TERMS; i++){
			data.poly[i] = coeff(polyzz,i);
		}

	cout << "(" << numgen << "," << numimp << ")";

	return data;
	
}

bool Vault::operator()(const folded& f1, const folded& f2) {
	        return method->compFold(f1,f2);
	    }