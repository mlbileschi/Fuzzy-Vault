#include "../inc/unitTests.h"
#include "../inc/ProbDistroGA.h"
#include "../inc/evaluatePerformance.h"
#include <iostream>
#include <stdexcept>
#include <vector>


bool Test0(){

	VaultMethod* method;
//	method = new VaultBF();
	method = new VaultTriangles();
	Vault* vault = new Vault(method);

	//Generate random secret
	vec_ZZ_p secret;
	secret.SetLength(POLYNOMIAL_TERMS);
	for(int i=0; i<POLYNOMIAL_TERMS; i++){
		secret[i] = random_ZZ_p();
	}

	vector<Minutiae> fingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/FVC2002/DB2/features/27_3.fp");
	for(int i=0; i<fingerprint.size(); i++){
		cout << fingerprint[i].m_nY << " ";
	}

	vault->lock(fingerprint, secret);
	polyResults result = vault->unlock(fingerprint);
	cout << " -" << result.score << "- ";

	if(result.score > 10.0){
		cout << "\nTest0: passed";
		return true;
	}
	cout << "\nTest0: failed";
	return false;
}

bool Test1(){
	if(1 == 0){
		cout << "\nTest1: passed";
		return true;
	}
	cout << "\nTest1: failed";
	return false;
}

bool evaluatePerformanceTest()
{
	cout<<"evaluatePerformanceTest:\t\t";
	cout.flush();
	vector<float> genuines;
	vector<float> impostors;
	genuines.push_back(-2);
	genuines.push_back(-1);
	impostors.push_back(0);
	impostors.push_back(1);
	bool toReturn = evaluatePerformance(impostors, genuines)==.5;

	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;
	return toReturn;
}

bool nextGenerationTest()
{
	cout<<"nextGenerationTest:\t\t\t";
	cout.flush();
	//arbitrary choices, except numFingerprints=100 (if we take a subset, the eer will change
	//from run to run. populationSize is low to speed it up
	int populationSize = 2;
	int genomeSize = 46;
	int numFingerprints = 100;
	float distroIncrement = .8;
	string database = "FVC2004/DB3";

	GenProps gp = nextGeneration(newGeneration(populationSize, genomeSize), populationSize, genomeSize, 
								numFingerprints, distroIncrement, database);

	bool toReturn = true;

	toReturn = toReturn && (gp.generation.size() == populationSize);
	for(int i = 0; i < gp.generation.size(); i++)
	{
		toReturn = toReturn && (gp.generation.at(i).size() == genomeSize);
	}

	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}

//FVC2006/DB1 is too low resolution to consider;
//there are readings with 0 minutiae
bool readAllFingerprintsTest()
{
	cout<<"readAllFingerprintsTest:\t\t";
	cout.flush();
	bool toReturn = true;

	//try good DB's
	for(int dbYear = 2002; dbYear<=2006; dbYear+=2)
	{
		for(int dbNum = 1; dbNum<=4; dbNum++)
		{
			//FVC2006/DB1 is too low resolution to consider;
			//there are readings with 0 minutiae
			if(dbYear == 2006 && dbNum ==1) continue;

			try
			{
				string goodDatabase = "FVC" + toString(dbYear) + "/DB" + toString(dbNum);
				vector<vector<vector<Minutiae> > > goodFingerprints = readAllFingerprints(goodDatabase);
				toReturn = toReturn && (goodFingerprints.size()==100);
				for(int i = 0; i<goodFingerprints.size(); i++)
				{
					toReturn = toReturn && (goodFingerprints.at(i).size()==8);
					for(int j = 0; j<8; j++)
					{
						toReturn = toReturn && (goodFingerprints.at(i).at(j).size()>0);
					}
				}
			}
			catch (out_of_range& oor)
			{
				toReturn = false;
			}
		}
	}
	string badDatabase = "kjasdf";

	try //bad DB name
	{
		vector<vector<vector<Minutiae> > > badFingerprints = readAllFingerprints(badDatabase);
		toReturn = toReturn && (badFingerprints.size()==100);
		for(int i = 0; i<badFingerprints.size(); i++)
		{
			toReturn = toReturn && (badFingerprints.at(i).size()==8);
			for(int j = 0; j<8; j++)
			{
				toReturn = toReturn && (badFingerprints.at(i).at(j).size()==0);
			}
		}
	}
	catch (out_of_range& oor)
	{
		toReturn = false;
	}


	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}

//makes sure the extern variable allFingerprints
//is read correctly, and other files have access to it
bool allFingerprintsReadCorrectlyTest()
{
	cout<<"allFingerprintsReadCorrectlyTest:\t";
	cout.flush();
	bool toReturn = true;

	toReturn = toReturn && (allFingerprints.size()==100);
	for(int i = 0; i<allFingerprints.size(); i++)
	{
		toReturn = toReturn && (allFingerprints.at(i).size()==8);
		for(int j = 0; j<8; j++)
		{
			toReturn = toReturn && (allFingerprints.at(i).at(j).size()>0);
		}
	}

	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}


//makes sure the extern variable allVaults
//is read correctly, and other files have access to it
//(only checks for the current DB)
bool allVaultsReadCorrectlyTest()
{
	cout<<"allVaultsReadCorrectlyTest:\t\t\t";
	cout.flush();
	bool toReturn = true;

	toReturn = toReturn && (allVaults.size()==100);
	for(int i = 0; i<allVaults.size(); i++)
	{
		toReturn = toReturn && (allVaults.at(i).size()==8);
	}
	for(int i = 0; i<allVaults.size(); i++)
	{
		for(int j = 0; j<allVaults.at(i).size(); j++)
		{
			toReturn = toReturn && ((allVaults.at(i).at(j)->unlock(allFingerprints.at(i).at(j))).score>5);
			//cout<<" "<<(allVaults.at(i).at(j)->unlock(allFingerprints.at(i).at(j))).score;
		}
	}

	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}

bool runDatabaseTest()
{
	cout<<"runDatabaseTest:\t\t\t";
	cout.flush();
	bool toReturn = true;

	int genomeSize = 20;
	string database = "FVC2004/DB4";
	int numFingerprints = 50;
	float distroIncrement = 1.1;
	vector<float> distro = newDistro(genomeSize);
	
	float eer = runDatabase(database, &distro, numFingerprints, distroIncrement);

	toReturn = toReturn && (eer < 1) && (eer>0);
	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}

void runUnitTests(){
	int count = 0;
	int pass = 0;

	count++;
	pass += evaluatePerformanceTest();
	count++;
	pass += nextGenerationTest();
	count++;
	pass += readAllFingerprintsTest();
	count++;
	pass += allFingerprintsReadCorrectlyTest();
	count++;
	pass += allVaultsReadCorrectlyTest();
	count++;
	pass += runDatabaseTest();


	cout << "\nResults: " << pass << " of " << count << " test(s) passed\n";
	if(count==pass){cout << "PASSED!!"<<endl;} else { cout << "Your code is bad and you should feel bad"<<endl;}
}
