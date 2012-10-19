#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "../inc/RunTest.h"
#include "../inc/EvaluatePerformance.h"
#include "../inc/Vault.h"
#include "../inc/IOUtils.h"
#include "../../matching_tjea/inc/matching.h"
#include "../../matching_tjea/inc/global_parameters.h"

using namespace std;

float root[5000000];
int zz1;
int zz2;
vector<vector<vector<folded> > > allFingerprintsFolded;
vector<vector<vector<vector<folded> > > > allFingerprintsFoldedTransformed;
vector<vector<vector<vector<folded> > > > allFingerprintsFoldedRotated;
vector<vector<vector<vector<minutia> > > > allFingerprintsRotatedMinu;
vector<vector<vector<minutia> > > allFingerprints;
 
//this method is here, and uses allFingerprints
//this is dummy data. the real DB will be read in in main
vector<vector<Vault *> > allVaults;// = dummy();

vector<vector<Vault * > > dummy()
{
	vector<vector<Vault *> >  allFP;
	return allFP;
}

void makeRoot(){
    for(int i=0; i<5000000; i++){
        //root.push_back(sqrt(to_float(i)));
        root[i] =sqrt(to_float(i));

       
    }
}

void makeAllVaults(string methodType, float distroIncrement)
{

	VaultMethod* method;

	if(methodType == "triangles")
		method = new VaultTriangles();
	else if(methodType == "bf")
		method = new VaultBF();
	else if(methodType == "new")
		method = new VaultLetsTrySomething();
	else //use the fast one
	{
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no match for method type. using default triangles."<<endl;
		method = new VaultTriangles();
	}	

	vector<vector<Vault *> > toReturn;
	for(int printNumber = 0; printNumber < allFingerprints.size(); printNumber++)
	{
		vector<Vault *> vaultsForPrintNumber;
		static int p=0;
		//cout << " " << ++p << " ";
		for(int readingNumber = 0; readingNumber<allFingerprints.at(printNumber).size(); readingNumber++)
		{
			//cout << vaultsForPrintNumber.size() << " ";
			Vault * vault = new Vault(method);
			//cout << vault << endl;
			
			vault->setDistroIncrement(distroIncrement);
			//Generate random secret
			vec_ZZ_p secret;
			secret.SetLength(POLYNOMIAL_TERMS);
			for(int i=0; i<POLYNOMIAL_TERMS; i++){
				secret[i] = random_ZZ_p();
			}

			vault->lock(allFingerprints.at(printNumber).at(readingNumber), secret);
			vaultsForPrintNumber.push_back(vault);
		}
		toReturn.push_back(vaultsForPrintNumber);
	}
	allVaults = toReturn;
}


void makeAllFolded(string methodType){

	//initialize the NTL stuff
	ZZ_p::init(FEILD_SIZE);

	VaultMethod* method;

	if(methodType == "triangles"){
		method = new VaultTriangles();}
	else if(methodType == "new"){ cout << "NEW" << flush;
		method = new VaultLetsTrySomething();
		return;}
	else if(methodType == "bf"){
		method = new VaultBF();}
	else //use the fast one
	{
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no match for method type. using default triangles."<<endl;
		method = new VaultTriangles();
	}	

	vector<vector<vector<folded> > > allFP;
	for(int fingerprintNumber=0; fingerprintNumber<allFingerprints.size(); fingerprintNumber++)
	{
		vector<vector<folded> > toPush;
		for(int readingNumber=0; readingNumber<allFingerprints[fingerprintNumber].size(); readingNumber++)
		{
			vector<folded> readFingerprint;
			folded point;
/*
			std::set<int> cc = method->minutiae2int(allFingerprints[fingerprintNumber][readingNumber]);
			for(std::set<int>::iterator iter=cc.begin(); iter != cc.end(); iter++){
				point.zint = *iter;
				readFingerprint.push_back(point);
			}
*/

			std::set<ZZ> cc = method->minutiae2ZZ(allFingerprints[fingerprintNumber][readingNumber]);
			for(std::set<ZZ>::iterator iter=cc.begin(); iter != cc.end(); iter++){
				point.z1 = to_ZZ_p(*iter);
				readFingerprint.push_back(point);
			}


			//cout << readFingerprint.size() << " ";

			toPush.push_back(readFingerprint);
		}
		allFP.push_back(toPush);
	}

	allFingerprintsFolded = allFP;
}



void makeAllFoldedTransformed(string methodType){

	VaultMethod* method;

	if(methodType == "triangles")
		//method = new VaultTriangles();
		return;
	else if(methodType == "bf")
		method = new VaultBF();
	else //use the fast one
	{
		return;
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no match for method type. using default triangles."<<endl;
		method = new VaultTriangles();
	}	


	int center_x = FP_IMAGE_WIDTH>>1;
	int center_y = FP_IMAGE_HEIGHT>>1;

	int theta_min, theta_max, x_min, x_max, y_min, y_max;

	theta_max = THETA_STEPS >> 1;
	x_max = X_STEPS >> 1;
	y_max = Y_STEPS >> 1;

	theta_min = THETA_STEPS & 1 ? -(THETA_STEPS>>1) : -(THETA_STEPS>>1) +1;
	x_min = X_STEPS & 1 ? -(X_STEPS>>1) : -(X_STEPS>>1) +1;
	y_min = Y_STEPS & 1 ? -(Y_STEPS>>1) : -(Y_STEPS>>1) +1;


	int d_x, d_y;
	float d_theta, rot_x, rot_y;
	float thisCos, thisSin;
	int n;

	vector<vector<vector<vector<folded> > > > allFP;

	for(int fingerprintNumber=0; fingerprintNumber<allFingerprints.size(); fingerprintNumber++)
	{
		vector<vector<vector<folded> > > allReadings;
		for(int readingNumber=0; readingNumber<allFingerprints[fingerprintNumber].size(); readingNumber++)
		{
			vector<vector<folded> > allTransformations;
			vector<minutia>* reading;
			reading = &allFingerprints[fingerprintNumber][readingNumber]; //alias to shorten the name
			//various rotations and translations
			n = reading->size();

			vector<minutia> testerRotated = vector<minutia>(n);
			vector<minutia> tester = vector<minutia>(n);
			minutia minut;
			vector<minutia>::iterator iter;
			vector<minutia>::iterator rotIter;

			for(short i=theta_min; i<=theta_max; i++){ 
				d_theta = i*THETA_INCREMENT;
				thisCos = cos(d_theta/57.3); //convert to radians (roughly)
				thisSin = sin(d_theta/57.3);

				for(int q=0; q<n; ++q){
					rot_x = center_x + thisCos*((*reading)[q].x-center_x) + thisSin*((*reading)[q].y-center_y);
					rot_y = center_y - thisSin*((*reading)[q].x-center_x) + thisCos*((*reading)[q].y-center_y);
					minut.x = (int)rot_x;
					minut.y = (int)rot_y;
					minut.theta = (*reading)[q].theta + (int)d_theta;
					testerRotated[q] = minut;
					tester[q] = minut;
				}	
	
				for(short j=x_min; j<=x_max; j++){

					d_x = (int)(j*X_INCREMENT);

					for(short k=y_min; k<=y_max; k++){

						d_y = (int)(k*Y_INCREMENT);
				
						iter = tester.begin();
				
						for(rotIter = testerRotated.begin(); rotIter<testerRotated.end() && iter<tester.end(); ++rotIter, ++iter){
							(*iter).x = (*rotIter).x + d_x;
							(*iter).y = (*rotIter).y + d_y; 
						}
						
						vector<folded> readFingerprint;
						std::set<ZZ> cc = method->minutiae2ZZ(tester);
						folded point;

						for(std::set<ZZ>::iterator iterz=cc.begin(); iterz != cc.end(); iterz++){
							point.z1 = to_ZZ_p(*iterz);
							readFingerprint.push_back(point);
						}

						allTransformations.push_back(readFingerprint);

					}
				}
				testerRotated = vector<minutia>(n);
				tester = vector<minutia>(n);
			}
			allReadings.push_back(allTransformations);
		}
		allFP.push_back(allReadings);
	}
	allFingerprintsFoldedTransformed = allFP;

//cout << endl << "prints: " << allFingerprintsFoldedTransformed.size();
//cout << endl << "readings: " << allFingerprintsFoldedTransformed[0].size();
cout << endl << "transformations: " << allFingerprintsFoldedTransformed[0][0].size();
//cout << endl << "minutiae: " << allFingerprintsFoldedTransformed[0][0][0].size();
//cout << endl << endl << flush;
cout << endl;
}






void makeAllFoldedRotated(string methodType){

	VaultMethod* method;

	if(methodType == "triangles" || methodType == "new"){
		//method = new VaultTriangles(); 
		return;
	}
	else if(methodType == "bf"){
		method = new VaultBF();
	}
	else //use the fast one
	{
		return;
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no match for method type. using default triangles."<<endl;
		method = new VaultTriangles();
	}	


	int center_x = FP_IMAGE_WIDTH>>1;
	int center_y = FP_IMAGE_HEIGHT>>1;

	int theta_min, theta_max, x_min, x_max, y_min, y_max;

	theta_max = THETA_STEPS >> 1;
	x_max = X_STEPS >> 1;
	y_max = Y_STEPS >> 1;

	theta_min = THETA_STEPS & 1 ? -(THETA_STEPS>>1) : -(THETA_STEPS>>1) +1;
	x_min = X_STEPS & 1 ? -(X_STEPS>>1) : -(X_STEPS>>1) +1;
	y_min = Y_STEPS & 1 ? -(Y_STEPS>>1) : -(Y_STEPS>>1) +1;


	int d_x, d_y;
	float d_theta, rot_x, rot_y;
	float thisCos, thisSin;
	int n;

	vector<vector<vector<vector<folded> > > > allFP;
	vector<vector<vector<vector<minutia> > > > allFPMinu; // Minu

	for(int fingerprintNumber=0; fingerprintNumber<allFingerprints.size(); fingerprintNumber++)
	{
		vector<vector<vector<folded> > > allReadings;
		vector<vector<vector<minutia> > > allReadingsMinu; // Minu
		for(int readingNumber=0; readingNumber<allFingerprints[fingerprintNumber].size(); readingNumber++)
		{
			vector<vector<folded> > allTransformations;
			vector<vector<minutia> > allTransformationsMinu; // Minu
			vector<minutia>* reading;
			reading = &allFingerprints[fingerprintNumber][readingNumber]; //alias to shorten the name
			//various rotations and translations
			n = reading->size();

			vector<minutia> testerRotated = vector<minutia>(n);
			vector<minutia> tester = vector<minutia>(n);
			minutia minut;
			vector<minutia>::iterator iter;
			vector<minutia>::iterator rotIter;

			for(short i=theta_min; i<=theta_max; i++){ 
				d_theta = i*THETA_INCREMENT;
				thisCos = cos(d_theta/57.3); //convert to radians (roughly)
				thisSin = sin(d_theta/57.3);

				vector<minutia> readFingerprintMinu; // Minu

				for(int q=0; q<n; ++q){
					rot_x = center_x + thisCos*((*reading)[q].x-center_x) - thisSin*((*reading)[q].y-center_y);
					rot_y = center_y + thisSin*((*reading)[q].x-center_x) + thisCos*((*reading)[q].y-center_y);
					minut.x = floor(rot_x+ .5);
					minut.y = floor(rot_y+ .5);
					minut.theta = floor((*reading)[q].theta + d_theta +.5);
					testerRotated[q] = minut;
					//tester[q] = minut;
					readFingerprintMinu.push_back(minut); // Minu
				}	
				
	
				/*for(short j=x_min; j<=x_max; j++){

					d_x = (int)(j*X_INCREMENT);

					for(short k=y_min; k<=y_max; k++){

						d_y = (int)(k*Y_INCREMENT);
				
						iter = tester.begin();
				
						for(rotIter = testerRotated.begin(); rotIter<testerRotated.end() && iter<tester.end(); ++rotIter, ++iter){
							(*iter).x = (*rotIter).x + d_x;
							(*iter).y = (*rotIter).y + d_y; 
						}
				*/
						vector<folded> readFingerprint;
						std::set<ZZ> cc = method->minutiae2ZZ(testerRotated);
						folded point;

						for(std::set<ZZ>::iterator iterz=cc.begin(); iterz != cc.end(); iterz++){
							point.z1 = to_ZZ_p(*iterz);
							readFingerprint.push_back(point);
						}

						allTransformations.push_back(readFingerprint);
						allTransformationsMinu.push_back(readFingerprintMinu); // Minu

					//}
				//}
				testerRotated = vector<minutia>(n);
				tester = vector<minutia>(n);
			}
			allReadings.push_back(allTransformations);
			allReadingsMinu.push_back(allTransformationsMinu); // Minu
		}
		allFP.push_back(allReadings);
		allFPMinu.push_back(allReadingsMinu); // Minu
	}
	allFingerprintsFoldedRotated = allFP;
	allFingerprintsRotatedMinu = allFPMinu; // Minu

//cout << endl << "prints: " << allFingerprintsFoldedTransformed.size();
//cout << endl << "readings: " << allFingerprintsFoldedTransformed[0].size();
//cout << endl << "transformations: " << allFingerprintsFoldedTransformed[0][0].size();
cout << endl << "Rotations: " << allFingerprintsFoldedRotated[0][0].size();
//cout << endl << "minutiae: " << allFingerprintsFoldedTransformed[0][0][0].size();
//cout << endl << endl << flush;
cout << endl;
}







vector<float> testGenuines(int fingerprintIndex, int readingNumber, Vault* vault)
{
	vector<float> toReturn;
	for(int i=readingNumber+1; i<=8; i++)
	{	
		
		if(vault->getMethodType()=="triangles")
			toReturn.push_back(vault->unlock(allFingerprintsFolded.at(fingerprintIndex-1).at(i-1)).score);	
		else if(vault->getMethodType()=="new")
			toReturn.push_back(vault->unlock(allFingerprints.at(fingerprintIndex-1).at(i-1)).score);
		else if (vault->getMethodType()=="bf"){
	//cout << " /" << flush;
			//toReturn.push_back(vault->unlock(allFingerprintsFoldedTransformed[fingerprintIndex-1][i-1]).score);
			//toReturn.push_back(vault->unlockRot(allFingerprintsFoldedRotated[fingerprintIndex-1][i-1]).score);
			//toReturn.push_back(vault->unlockRot(allFingerprintsRotatedMinu[fingerprintIndex-1][i-1]).score);
			vector<pair<int,int> > matchedPairs;
//toReturn.push_back(match_tjea(allFingerprints[fingerprintIndex-1][readingNumber-1], allFingerprints[fingerprintIndex-1][i-1], matchedPairs));
//toReturn.push_back(match_tjea(vault->vault, allFingerprints[fingerprintIndex-1][i-1], matchedPairs));
			

			match_tjea(allFingerprints[fingerprintIndex-1][readingNumber-1], allFingerprints[fingerprintIndex-1][i-1], matchedPairs);
			vector<minutia> testPrint = 
			 translate(allFingerprints[fingerprintIndex-1][readingNumber-1], allFingerprints[fingerprintIndex-1][i-1], matchedPairs);
			toReturn.push_back(vault->unlock(testPrint).score);


			//toReturn.push_back(vault->unlock(allFingerprints.at(fingerprintIndex-1).at(i-1)).score);
		}
	}
	return toReturn;
}


vector<float> testImpostors(int fingerprintIndex, Vault* vault)
{
	vector<float> toReturn;
	for(int i=0; i<100; i++)
	{
		if(i==fingerprintIndex-1) 
			continue;
		if(vault->getMethodType()=="triangles")
			toReturn.push_back(vault->unlock(allFingerprintsFolded.at(i).at(0)).score);
		else if(vault->getMethodType()=="new")
			toReturn.push_back(vault->unlock(allFingerprints.at(i).at(0)).score);
		else if (vault->getMethodType()=="bf"){
	//cout << " ." << flush;
			//toReturn.push_back(vault->unlock(allFingerprintsFoldedTransformed[i][0]).score);
			//toReturn.push_back(vault->unlockRot(allFingerprintsFoldedRotated[i][0]).score);
			//toReturn.push_back(vault->unlockRot(allFingerprintsRotatedMinu[i][0]).score);
			//toReturn.push_back(-5);
			vector<pair<int,int> > matchedPairs;
			//toReturn.push_back((float)match_tjea(allFingerprints[fingerprintIndex-1][0], allFingerprints[i][0], matchedPairs));
			//toReturn.push_back((float)match_tjea(vault->vault, allFingerprints[i][0], matchedPairs));


			match_tjea(allFingerprints[fingerprintIndex-1][0], allFingerprints[i][0], matchedPairs);
			vector<minutia> testPrint = 
			 translate(allFingerprints[fingerprintIndex-1][0], allFingerprints[i][0], matchedPairs);
			toReturn.push_back(vault->unlock(testPrint).score);


			//toReturn.push_back(vault->unlock(allFingerprints[i][0]).score);
			//toReturn.push_back(vault->unlock(allFingerprints.at(i).at(0)).score);
		}
	}
	return toReturn;
}

//@param database is like "FVC2002/DB1"
//@param numFingerprints needs to be >=1, <=100 when given
//default value is 100
float runDatabase(string database, vector<float>* distro, int numFingerprints, float distroIncrement)
{
	vector<float> genuineScores;
	vector<float> impostorScores;

	int randStart = rand()%100; //rand number between 0 and 99
	int fingerprintIndex = 0;
	int numExamined = 0;


	//if you plan on running this on different databases than the FVC200[246],
	//you're going to need to change this 100, and the 8's
	while(numExamined<numFingerprints)
	{

		fingerprintIndex = ((randStart + numExamined) % 100)+1;
		//cout << endl << "fingerprint Index: " << fingerprintIndex << endl;
		for(int readingNumber = 1; readingNumber<=8; readingNumber++)
		{
			Vault * vault = allVaults.at(fingerprintIndex-1).at(readingNumber-1);
			vault->setDistro(*distro);
			vector<float> newGenScores = testGenuines(fingerprintIndex, readingNumber, vault);
			genuineScores.insert(genuineScores.end(), newGenScores.begin(), newGenScores.end() );

			//only test the first readings for impostors.
			//we have a lot more impostors than genuines, 
			//so we must choose which ones to test
			if(readingNumber==1)
			{
				//cout << endl;
				vector<float> newImpScores = testImpostors(fingerprintIndex, vault);
				impostorScores.insert(impostorScores.end(), newImpScores.begin(), newImpScores.end());
				//cout << endl;
			}

		}
	
		numExamined++;
	}

	return evaluatePerformance(impostorScores, genuineScores);
}


double match_tjea(vector<minutia> finger1, vector<minutia> finger2, vector<pair<int,int> > &matchingPairs)
{

	FPTemplate f1;
	FPTemplate f2;

	ptrFPTemplate fp1;
	fp1 = &f1;
	ptrFPTemplate fp2;
	fp2 = &f2;

	fp1->nCnt = finger1.size();
	for(int i=0; i<fp1->nCnt; i++){
		fp1->marr[i].m_nX = finger1[i].x;
		fp1->marr[i].m_nY = finger1[i].y;
		fp1->marr[i].m_nTheta = finger1[i].theta;
	}

	fp2->nCnt = finger2.size();
	for(int i=0; i<fp2->nCnt; i++){
		fp2->marr[i].m_nX = finger2[i].x;
		fp2->marr[i].m_nY = finger2[i].y;
		fp2->marr[i].m_nTheta = finger2[i].theta;
	}


	MatchResult result;
	MatchResultEx unused;

	CUBS_MatchFeaturesInternal(fp1, fp2, &result, &unused, matchingPairs);

//	cout << endl << endl << "!!---------------    " << result.similarity << "   -----------------!!" << endl << endl;

	return result.similarity;

}




// rotates and translates test onto reference so the matching points line up well

vector<minutia> translate(const vector<minutia> &reference, const vector<minutia> &test, vector<pair<int, int> > &matchingPoints)
{
	vector<minutia> toReturn;
	minutia minut;

	if(matchingPoints.size() < 1)
	{	
		minutia minut;
		for(int i=0; i<test.size(); i++){
			minut.x = test[i].x;
			minut.y = test[i].y;
			minut.theta = test[i].theta;
			toReturn.push_back(minut);
		}
		return toReturn;
	}	
	float rotation;
	int xShift;
	int yShift;
	
	minutia ref = reference[matchingPoints[0].second];
	minutia r1;// = reference[matchingPoints[1].second];

	minutia t0 = test[matchingPoints[0].first];
	minutia t1;

	xShift = ref.x - t0.x;
	yShift = ref.y - t0.y;
	

	for(int i=0; i<matchingPoints.size(); i++)
	{
		r1 = reference[matchingPoints[i].second];
		t1 =  test[matchingPoints[i].first];

		int dtheta = r1.theta - t1.theta;
		if(dtheta > 180){dtheta -= 360;}
		if(dtheta < -180){dtheta += 360;}

		//cout << " " << dtheta << " ";
	}


	for(int i=0; i<test.size(); i++){
		minut.x = test[i].x + xShift;
		minut.y = test[i].y + yShift;
		minut.theta = test[i].theta;		
		toReturn.push_back(minut);
	}

	t0 = toReturn[matchingPoints[0].first];
	t1 = toReturn[matchingPoints[1].first];

	int center_x = ref.x;
	int center_y = ref.y;

	double sum = 0.0;
	double px,py,radiansp,qx,qy,radiansq,angle,d_rad;
//cout << " |-| ";
	for(int i=1; i<matchingPoints.size(); i++)
	{
		r1 = reference[matchingPoints[i].second];
		t1 =  toReturn[matchingPoints[i].first];

		px = r1.x - ref.x;
		py = r1.y - ref.y;
		radiansp = atan2(py,px);

		qx = t1.x - ref.x;
		qy = t1.y - ref.y;
		radiansq = atan2(qy,qx);
	//cout << " " << radiansp - radiansq;
		d_rad = radiansp - radiansq;
		if(d_rad > 3.1415){d_rad -= 2.0*3.1415;}
		if(d_rad < -3.1415){d_rad += 2.0*3.1415;}
		sum += d_rad;
	}

	angle = sum/(matchingPoints.size()-1);

//cout << "angle: ";
//cout << angle << " ";

	double thisCos = cos(angle); 
	double thisSin = sin(angle);

	//rotate em
	for(int i=1; i<test.size(); i++){
		toReturn[i].x = floor(center_x + thisCos*(toReturn[i].x-center_x) - thisSin*(toReturn[i].y-center_y) +.5);
		toReturn[i].y = floor(center_y + thisSin*(toReturn[i].x-center_x) + thisCos*(toReturn[i].y-center_y) +.5);
		toReturn[i].theta += floor(angle*(180/3.1415) +.5);
	}

	for(int i=0; i<matchingPoints.size(); i++){
		r1 = reference[matchingPoints[i].second];
		t1 =  toReturn[matchingPoints[i].first];

		float Root = sqrt(pow(r1.x-t1.x,2.0) + pow(r1.y-t1.y,2.0));
		//cout << "!!! " << Root << " !!!";
	}

	int shiftx,shifty;
	int sumx=0;
	int sumy=0;

	for(int i=0; i<matchingPoints.size(); i++)
	{
		r1 = reference[matchingPoints[i].second];
		t1 =  toReturn[matchingPoints[i].first];

		sumx += r1.x - t1.x;
		sumy += r1.y - t1.y;
	//cout << "()()() " << sumy << " ";
	}
	int n = matchingPoints.size();
	//cout << matchingPoints.size() << "  ";
	shiftx = sumx/n;
	shifty = sumy/n;

	//cout << "^^^^^^^^^ " << shiftx << " " << shifty << "  ";

	for(int i=1; i<test.size(); i++){
		toReturn[i].x += shiftx;
		toReturn[i].y += shifty;
	}


	for(int i=0; i<matchingPoints.size(); i++){
		r1 = reference[matchingPoints[i].second];
		t1 =  toReturn[matchingPoints[i].first];

		float Root = sqrt(pow(r1.x-t1.x,2.0) + pow(r1.y-t1.y,2.0));
		//cout << " " << Root << " ";
	}



	return toReturn;
}








