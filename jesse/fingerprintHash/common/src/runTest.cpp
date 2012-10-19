#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "../inc/runTest.h"
#include "../inc/evaluatePerformance.h"
#include "../inc/vault.h"
#include "../inc/IOUtils.h"

using namespace std;

//this method is from ioutils
//this is dummy data. the real DB will be read in in main
vector<vector<vector<folded> > > allFingerprintsFolded;
vector<vector<vector<vector<folded> > > > allFingerprintsFoldedTransformed;
vector<vector<vector<Minutiae> > > allFingerprints;// = readAllFingerprints("FVC2002/DB1"); 
//this method is here, and uses allFingerprints
//they have been locked
//this is dummy data. the real DB will be read in in main
vector<vector<Vault *> > allVaults;// = dummy();

vector<vector<Vault * > > dummy()
{
	vector<vector<Vault *> >  allFP;
	return allFP;
}
//needs to be made configurable by the actual method you want to choose
void makeAllVaults(string methodType, float distroIncrement)
{


	VaultMethod* method;

	if(methodType == "triangles")
		method = new VaultTriangles();
	else if(methodType == "bf")
		method = new VaultBF();
	else //use the fast one
	{
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no method specified. using default triangles"<<endl;
		method = new VaultTriangles();
	}	

	vector<vector<Vault *> > toReturn;
	for(int printNumber = 0; printNumber < allFingerprints.size(); printNumber++)
	{
		vector<Vault *> vaultsForPrintNumber;
		for(int readingNumber = 0; readingNumber<allFingerprints.at(printNumber).size(); readingNumber++)
		{
			Vault * vault = new Vault(method);
			vault->setDistroIncrement(distroIncrement);
			//Generate random secret
			vec_ZZ_p secret;
			secret.SetLength(POLYNOMIAL_TERMS);
			for(int i=0; i<POLYNOMIAL_TERMS; i++){
				secret[i] = random_ZZ_p();
			}

			vault->lock(allFingerprints.at(printNumber).at(readingNumber), secret);
			vaultsForPrintNumber.push_back(vault);
			//delete vault; //causes segfault if included, fyi
		}
		toReturn.push_back(vaultsForPrintNumber);
	}
//	delete method; //causes segfault if included, fyi
	allVaults = toReturn;
}


void makeAllFolded(string methodType){

	//initialize the NTL stuff
	ZZ_p::init(FEILD_SIZE);

	VaultMethod* method;

	if(methodType == "triangles")
		method = new VaultTriangles();
	else if(methodType == "bf")
		method = new VaultBF();
	else //use the fast one
	{
		cout<<"METHOD TYPE == |"<<methodType<<"|"<<endl;
		cout<<"no method specified. using default triangles"<<endl;
		method = new VaultTriangles();
	}	

	vector<vector<vector<folded> > > allFP;
	for(int fingerprintNumber=0; fingerprintNumber<allFingerprints.size(); fingerprintNumber++)
	{
		vector<vector<folded> > toPush;
		for(int readingNumber=0; readingNumber<allFingerprints[fingerprintNumber].size(); readingNumber++)
		{
			vector<folded> readFingerprint;
			std::set<ZZ> cc = method->minutiae2ZZ(allFingerprints[fingerprintNumber][readingNumber]);
			folded point;

			for(std::set<ZZ>::iterator iter=cc.begin(); iter != cc.end(); iter++){
				point.z1 = to_ZZ_p(*iter);
				readFingerprint.push_back(point);
			}

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
		cout<<"no method specified. using default triangles"<<endl;
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
			vector<Minutiae>* reading;
			reading = &allFingerprints[fingerprintNumber][readingNumber]; //alias to shorten the name
			//various rotations and translations
			n = reading->size();

			vector<Minutiae> testerRotated = vector<Minutiae>(n);
			vector<Minutiae> tester = vector<Minutiae>(n);
			Minutiae minut;
			vector<Minutiae>::iterator iter;
			vector<Minutiae>::iterator rotIter;

			for(short i=theta_min; i<=theta_max; i++){ 
				d_theta = i*THETA_INCREMENT;
				thisCos = cos(d_theta/57.3); //convert to radians (roughly)
				thisSin = sin(d_theta/57.3);

				for(int q=0; q<n; ++q){
					rot_x = center_x + thisCos*((*reading)[q].m_nX-center_x) + thisSin*((*reading)[q].m_nY-center_y);
					rot_y = center_y - thisSin*((*reading)[q].m_nX-center_x) + thisCos*((*reading)[q].m_nY-center_y);
					minut.m_nX = (int)rot_x;
					minut.m_nY = (int)rot_y;
					minut.m_nTheta = (*reading)[q].m_nTheta + (int)d_theta;
					testerRotated[q] = minut;
					tester[q] = minut;
				}	
	
				for(short j=x_min; j<=x_max; j++){

					d_x = (int)j*X_INCREMENT;

					for(short k=y_min; k<=y_max; k++){

						d_y = (int)k*Y_INCREMENT;
				
						iter = tester.begin();
				
						for(rotIter = testerRotated.begin(); rotIter<testerRotated.end() && iter<tester.end(); ++rotIter, ++iter){
							(*iter).m_nX = (*rotIter).m_nX + d_x;
							(*iter).m_nY = (*rotIter).m_nY + d_y; 
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
				testerRotated = vector<Minutiae>(n);
				tester = vector<Minutiae>(n);
			}
			allReadings.push_back(allTransformations);
		}
		allFP.push_back(allReadings);
	}
	allFingerprintsFoldedTransformed = allFP;

//cout << "\nprints: " << allFingerprintsFoldedTransformed.size();
//cout << "\nreadings: " << allFingerprintsFoldedTransformed[0].size();
cout << endl << "transformations: " << allFingerprintsFoldedTransformed[0][0].size();
//cout << "\nminutiae: " << allFingerprintsFoldedTransformed[0][0][0].size();
cout << endl << endl << flush;
}


vector<float> testGenuines(int fingerprintIndex, int readingNumber, Vault* vault)
{
	vector<float> toReturn;
	for(int i=readingNumber+1; i<=8; i++)
	{	
		//toReturn.push_back(vault->unlock(allFingerprintsFolded.at(fingerprintIndex-1).at(i-1)).score);
		//toReturn.push_back(vault->unlock(allFingerprints.at(fingerprintIndex-1).at(i-1)).score);
		toReturn.push_back(vault->unlock(allFingerprintsFoldedTransformed[fingerprintIndex-1][i-1]).score);
	}
	return toReturn;
}

vector<float> testImpostors(int fingerprintIndex, Vault* vault)
{
	vector<float> toReturn;
	for(int i=0; i<100; i++)
	{
		if(i==fingerprintIndex) continue;
		//toReturn.push_back(vault->unlock(allFingerprintsFolded.at(i).at(0)).score);
		//toReturn.push_back(vault->unlock(allFingerprints.at(i).at(0)).score);
		toReturn.push_back(vault->unlock(allFingerprintsFoldedTransformed[i][0]).score);
	}
	return toReturn;
}

//@param database is like "FVC2002/DB1"
//@param numFingerprints is optional. needs to be >=1, <=100 when given
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
	while(numExamined<=numFingerprints)
	{
		
		fingerprintIndex = ((randStart + numExamined) % 100)+1;

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
				vector<float> newImpScores = testImpostors(fingerprintIndex, vault);
				impostorScores.insert(impostorScores.end(), newImpScores.begin(), newImpScores.end());
			}

		}
	
		numExamined++;
	}

	return evaluatePerformance(impostorScores, genuineScores);
	//return 0.0;
}
