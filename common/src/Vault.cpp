
#include "../inc/Vault.h"

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
	distro = vector<float>();
	//Vault::distroIncrement = 1.0;
	//distro = vector<float>();
	//Vault::distro.push_back(1.0);
	method->getDefaultDistro(&distro, &distroIncrement);
	computeMaxDistro();
}


string Vault::getMethodType(){
	return method->getMethodType();
}

//void Vault::lock(vector<folded> hRef, vec_ZZ_p secret){


void Vault::lock(vector<minutia> hRef, vec_ZZ_p secret){

	//Add genuine points
	std::set<ZZ> zValues = method->minutiae2ZZ(hRef);
	//std::set<int> zValuesint = method->minutiae2int(hRef);
	folded fold;
	fold.chaff = false;		


	//std::set<int>::iterator iterint=zValuesint.begin();
	for(std::set<ZZ>::iterator iter=zValues.begin(); iter != zValues.end() ; iter++){
		fold.z1 = to_ZZ_p(*iter);
		fold.f1 = function(to_ZZ_p(*iter), secret);
		vault.push_back(fold);
	}



	//std::set<int>::iterator iterint=zValuesint.begin();
	//for(std::set<int>::iterator iterint=zValuesint.begin(); iterint != zValuesint.end() ; iterint++){
	//	fold.zint = *iterint;
		//fold.f1 = to_ZZ_p(0);
		//vault.push_back(fold);
	//}

	//Add chaff points

	for(int i=0; i<NUM_CHAFF_POINTS; i++){

		folded chaff;
		chaff.chaff = true;

		chaff.z1 = random_ZZ_p();
		chaff.zint = rand();
		chaff.f1 = random_ZZ_p();

		//while((zValues.count(rep(chaff.z1)) == 1) || (chaff.f1 == function(chaff.z1, secret)) ){
		//	chaff.z1 = random_ZZ_p();
		//}

		//zValues.insert(rep(chaff.z1));
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
cout << endl;
for(int i=0; i<vault.size(); i++){
	cout << vault[i].chaff << " " << vault[i].z1 << " | ";
}
cout << endl << endl;
*/

/*
	for(int q=0; q<vault.size(); q++){
		cout << (vault)[q].z1 << endl;
	}
	cout << endl;
*/
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////

polyResults Vault::unlock(const vector<vector<folded> >& hTst){

	polyResults data;
	std::set<int> indexPoints;

	indexPoints = compareExtractBF(&data, hTst);
	keyInLock(&data, indexPoints);
	return data;
}
polyResults Vault::unlock(const vector<folded>& hTst){

	polyResults data;
	std::set<int> indexPoints;

	indexPoints = method->INV ? compareExtract(&data, hTst) : compareExtractBF(&data, hTst);
	keyInLock(&data, indexPoints);
	return data;
}
//I'm seeing triple!
polyResults Vault::unlock(const vector<minutia>& hTst){
//cout << " 1 " << hTst.size() << "  " << vault.size();
	polyResults data;
	std::set<int> indexPoints;

	indexPoints = method->INV ? compareExtract(&data, hTst) : compareExtractBF(&data, hTst);
	keyInLock(&data, indexPoints);
	return data;
}
polyResults Vault::unlockRot(vector<vector<folded> >& hTst){

	polyResults data;
	std::set<int> indexPoints;

	indexPoints = compareExtractBFRot(&data, hTst);
	keyInLock(&data, indexPoints);
	return data;
}
polyResults Vault::unlockRot(vector<vector<minutia> >& hTst){

	polyResults data;
	std::set<int> indexPoints;

	indexPoints = compareExtractBFRot(&data, hTst);
	keyInLock(&data, indexPoints);
	return data;
}

void Vault::keyInLock(polyResults* data, const std::set<int>& indexPoints){
//cout << " " << indexPoints.size() << " ";
	int numPoints = indexPoints.size();

	vec_ZZ_p xvec, fofxvec;
	xvec.SetLength(numPoints);
	fofxvec.SetLength(numPoints);

	int numgen=0, numimp=0;
	std::set<int>::iterator it;
	int i=0;
	for(it = indexPoints.begin(); it != indexPoints.end(); it++, i++ ) {
		xvec[i] = vault[*it].z1;
		fofxvec[i] = vault[*it].f1;
		vault[*it].chaff ? numimp++ : numgen++;
	}  

	ZZ_pX polyzz;

	data->poly.SetLength(POLYNOMIAL_TERMS);
	if(((numPoints-POLYNOMIAL_TERMS)>>1) > 0){
		//bool t = WelchBerlekamp (polyzz, xvec , fofxvec ,  (numPoints-POLYNOMIAL_TERMS)>>1);
	}

		for(int i=0; i<POLYNOMIAL_TERMS; i++){
		//	data->poly[i] = coeff(polyzz,i);
		}

	//cout << "(" << numgen << "," << numimp << ")" << data->score << flush;
	
}


bool Vault::operator()(const folded& f1, const folded& f2) {
	        return method->compFold(f1,f2);
	    }



std::set<int> Vault::compareExtract(polyResults* data, const vector<folded>& candidate){
//cout << " 4 ";
	data->score = 0.0;
	std::set<int> indexPoints;

	//int base = 0;
	int n = vault.size();
	int m = candidate.size();

	float dist;

	float minDist;
	int minIndex;
	float rando;
	int distroIndex;
//cout << "n " << n << "m " << m << " ";
	for(int i=0; i<m; i++){
		//minDist = maxDistro + method->getDistroIncrement(); // make it a little bigger than the max
		minDist = maxDistro;
		minIndex = -1;
		for(int k=0; k<n; k++){
//cout << " 5 ";
			dist = method->distance(candidate[i].z1, vault[k].z1);
			//dist = method->distance(candidate[i].zint, vault[k].zint);

			if(dist<=minDist){
				minDist = dist;
				minIndex = k;
			}
		}//end k
		if(minIndex != -1){
			//cout << minDist;
			rando = ((float)rand())/((float)RAND_MAX);
			distroIndex = (int)ceil(minDist/distroIncrement);
			if(!indexPoints.count(minIndex)){data->score += (vault[minIndex].chaff ? -distro[distroIndex] : distro[distroIndex]);}
			//data->score += vault[minIndex].chaff ? -1 : 1;
			if(rando < distro[distroIndex]){
				indexPoints.insert(minIndex);
			//	data->score += vault[minIndex].chaff ? -1 : 1;
			}
		}
	}//end i
	return indexPoints;
}



std::set<int> Vault::compareExtract(polyResults* data, const vector<minutia>& hTst){

	vector<folded> candidate;
	std::set<ZZ> cc = method->minutiae2ZZ(hTst);
	folded point;

	for(std::set<ZZ>::iterator iter=cc.begin(); iter != cc.end(); iter++){
		point.z1 = to_ZZ_p(*iter);
		candidate.push_back(point);
	}

	return compareExtract(data, candidate);
}



std::set<int> Vault::compareExtractint(polyResults* data, const vector<minutia>& hTst){
//cout << " 3 ";
	//cout << "# " << hTst.size() << " ";
	vector<folded> candidate;
	std::set<int> cc = method->minutiae2int(hTst);
	folded point;

	for(std::set<int>::iterator iter=cc.begin(); iter != cc.end(); iter++){
		point.zint = *iter;
		candidate.push_back(point);
	}

	return compareExtract(data, candidate);
}


//Doesn't do anything
std::set<int> Vault::compareExtractBF(polyResults* data, const vector<folded>& translations){
	std::set<int> sol;
	return sol;
}

std::set<int> Vault::compareExtractBF(polyResults* data, const vector<vector<folded> >& translations){
//cout << "." << flush;
	std::set<int> sol;
	polyResults tempResults;
	float maxi = -100000.0;
	int n = translations.size();

	for(int i=0; i<n; ++i){
		std::set<int> res = compareExtract(&tempResults, translations[i]);
		if(tempResults.score > maxi){
			maxi=tempResults.score;
			sol = res;
		}
	}

	data->score = maxi;

	return sol;
}

std::set<int> Vault::compareExtractBFRot(polyResults* data,  vector<vector<folded> >& rotations){
//cout << "." << flush;
	std::set<int> sol;
	polyResults tempResults;
	float maxi = -100000.0;
	int n = rotations.size();

	int x_min, x_max, y_min, y_max;

	x_max = X_STEPS >> 1;
	y_max = Y_STEPS >> 1;

	x_min = X_STEPS & 1 ? -(X_STEPS>>1) : -(X_STEPS>>1) +1;
	y_min = Y_STEPS & 1 ? -(Y_STEPS>>1) : -(Y_STEPS>>1) +1;
	

	vector<folded>::iterator iter;
	vector<folded>::iterator rotIter;

	int d_x, d_y;

	for(int i=0; i<n; ++i){

		vector<folded> toComp = vector<folded>(rotations[i].size());
		
		for(short j=x_min; j<=x_max; j++){

			d_x = (int)(j*X_INCREMENT);

			for(short k=y_min; k<=y_max; k++){

				d_y = (int)(k*Y_INCREMENT);

				iter = toComp.begin();
				
				for(rotIter = rotations[i].begin(); rotIter<rotations[i].end() && iter<toComp.end(); ++rotIter, ++iter){
					(*iter).z1 = (*rotIter).z1 + to_ZZ_p(d_x << 10) + to_ZZ_p(d_y << 4); //unsafe if >p
				}	
				
				std::set<int> res = compareExtract(&tempResults, toComp);
				if(tempResults.score > maxi){
					maxi=tempResults.score;
					sol = res;
				}

			}//end k

		}//end j

	}//end i

	data->score = maxi;

	return sol;
}


std::set<int> Vault::compareExtractBFRot(polyResults* data,  vector<vector<minutia> >& hTst){
//cout << "." << flush;
	std::set<int> sol;

	//int center_x = FP_IMAGE_WIDTH>>1;
	//int center_y = FP_IMAGE_HEIGHT>>1;

	float max = 0;
	polyResults tempResults;

	//Test various rotations and translations
	//int theta_min, theta_max;
	int x_min, x_max, y_min, y_max;

	//theta_max = THETA_STEPS >> 1;
	x_max = X_STEPS >> 1;
	y_max = Y_STEPS >> 1;

	//theta_min = THETA_STEPS & 1 ? -(THETA_STEPS>>1) : -(THETA_STEPS>>1) +1;
	x_min = X_STEPS & 1 ? -(X_STEPS>>1) : -(X_STEPS>>1) +1;
	y_min = Y_STEPS & 1 ? -(Y_STEPS>>1) : -(Y_STEPS>>1) +1;

	float maxi = -10000000.0;

	int d_x, d_y;
	//float d_theta, rot_x, rot_y;
	//float thisCos, thisSin;
	int n = hTst.size();

	vector<minutia> testerRotated = vector<minutia>(n);
	vector<minutia> tester = vector<minutia>(n);
	minutia minut;
	vector<minutia>::iterator iter;
	vector<minutia>::iterator rotIter;


	for(short i=0; i<n; i++){ 
//cout << endl<<"i: " << i;
		//d_theta = i*THETA_INCREMENT;
		//thisCos = cos(d_theta/57.3); //convert to radians (roughly)
		//thisSin = sin(d_theta/57.3);

		for(int q=0; q<n; ++q){
			//rot_x = center_x + thisCos*(hTst.at(q).x-center_x) + thisSin*(hTst.at(q).y-center_y);
			//rot_y = center_y - thisSin*(hTst.at(q).x-center_x) + thisCos*(hTst.at(q).y-center_y);
			//minut.x = (int)rot_x;
			//minut.y = (int)rot_y;
			//minut.theta = hTst.at(q).theta + (int)d_theta;
			testerRotated[q] = hTst[i][q];
			tester[q] = hTst[i][q];
		}	
	
		for(short j=x_min; j<=x_max; j++){
//cout << endl<<"j: " << j;

			d_x = (int)(j*X_INCREMENT);

			for(short k=y_min; k<=y_max; k++){
//cout << endl<<"k: " << k;
		
				d_y = (int)(k*Y_INCREMENT);
				
				iter = tester.begin();
				for(rotIter = testerRotated.begin(); rotIter<testerRotated.end() && iter<tester.end(); ++rotIter, ++iter){
					(*iter).x = (*rotIter).x + d_x;
					(*iter).y = (*rotIter).y + d_y; 
				}

				//std::set<int> res = compareExtract(&tempResults, tester);
				std::set<int> res = compareExtractint(&tempResults, tester);
				//tester = vector<minutia>(150);
				//tester.clear();
				//while(tester.size()){tester.pop_back();}
				if(tempResults.score > maxi){
					maxi=tempResults.score;
					sol = res;
				}
			}
		}
		testerRotated = vector<minutia>(n);
		tester = vector<minutia>(n);
	}

	data->score = maxi;
	return sol;
}


std::set<int> Vault::compareExtractBF(polyResults* data, const vector<minutia>& hTst){
//cout << "." << flush;
	std::set<int> sol;
//cout << " 2 ";
	int center_x = FP_IMAGE_WIDTH>>1;
	int center_y = FP_IMAGE_HEIGHT>>1;

	float max = 0;
	polyResults tempResults;

	//Test various rotations and translations
	int theta_min, theta_max, x_min, x_max, y_min, y_max;

	theta_max = THETA_STEPS >> 1;
	x_max = X_STEPS >> 1;
	y_max = Y_STEPS >> 1;

	theta_min = THETA_STEPS & 1 ? -(THETA_STEPS>>1) : -(THETA_STEPS>>1) +1;
	x_min = X_STEPS & 1 ? -(X_STEPS>>1) : -(X_STEPS>>1) +1;
	y_min = Y_STEPS & 1 ? -(Y_STEPS>>1) : -(Y_STEPS>>1) +1;

	float maxi = -10000000.0;

	int d_x, d_y;
	float d_theta, rot_x, rot_y;
	float thisCos, thisSin;
	int n = hTst.size();

	vector<minutia> testerRotated = vector<minutia>(n);
	vector<minutia> tester = vector<minutia>(n);
	minutia minut;
	vector<minutia>::iterator iter;
	vector<minutia>::iterator rotIter;

	for(short i=theta_min; i<=theta_max; i++){ 
//cout << endl<<"i: " << i;
		d_theta = i*THETA_INCREMENT;
		thisCos = cos(d_theta/57.29747); //convert to radians (roughly)
		thisSin = sin(d_theta/57.29747);

		for(int q=0; q<n; ++q){
			rot_x = center_x + thisCos*(hTst.at(q).x-center_x) + thisSin*(hTst.at(q).y-center_y);
			rot_y = center_y - thisSin*(hTst.at(q).x-center_x) + thisCos*(hTst.at(q).y-center_y);
			minut.x = (int)rot_x;
			minut.y = (int)rot_y;
			minut.theta = hTst.at(q).theta + (int)d_theta;
			testerRotated.at(q) = minut;
			tester.at(q) = minut;
		}	
	
		for(short j=x_min; j<=x_max; j++){
//cout << endl<<"j: " << j;

			d_x = (int)(j*X_INCREMENT);

			for(short k=y_min; k<=y_max; k++){
//cout << endl<<"k: " << k;
		
				d_y = (int)(k*Y_INCREMENT);
				
				iter = tester.begin();
				
				for(rotIter = testerRotated.begin(); rotIter<testerRotated.end() && iter<tester.end(); ++rotIter, ++iter){
					(*iter).x = (*rotIter).x + d_x;
					(*iter).y = (*rotIter).y + d_y; 
				}

				std::set<int> res = compareExtract(&tempResults, tester);
				//std::set<int> res = compareExtractint(&tempResults, tester);
				//tester = vector<minutia>(150);
				//tester.clear();
				//while(tester.size()){tester.pop_back();}
				if(tempResults.score > maxi){
					maxi=tempResults.score;
					sol = res;
				}
			}
		}
		testerRotated = vector<minutia>(n);
		tester = vector<minutia>(n);
	}

	data->score = maxi;
	return sol;
}



void Vault::setDistro(const vector<float>& d){distro = d; computeMaxDistro(); distro.push_back(0.0);}
std::vector<float> Vault::getDistro(){return distro;}

void Vault::setDistroIncrement(const float& inc){distroIncrement = inc; computeMaxDistro();}
float Vault::getDistroIncrement(){return distroIncrement;}

void Vault::computeMaxDistro(){

	maxDistro = distro.size();
	float j = 0;
	while(j==0){
		maxDistro--;
		j = distro[(int)maxDistro];
	}

	maxDistro *= distroIncrement;


}


