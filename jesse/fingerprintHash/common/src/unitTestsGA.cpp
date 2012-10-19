#include "../inc/unitTests.h"
#include "../inc/ProbDistroGA.h"
#include "../inc/evaluatePerformance.h"
#include <iostream>
#include <stdexcept>
#include <vector>


bool evaluatePerformanceTest()
{
	cout<<"evaluatePerformanceTest:\t\t\t";
	cout.flush();

	readAllFingerprints("FVC2002/DB1");
	makeAllVaults("triangles", 1.0);

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
	cout<<"nextGenerationTest:\t\t\t\t";
	cout.flush();

	//arbitrary choices, except numFingerprints=100 (if we take a subset, the eer will change
	//from run to run. populationSize is low to speed it up
	int populationSize = 2;
	int genomeSize = 46;
	int numFingerprints = 100;
	float distroIncrement = .8;
	string database = "FVC2004/DB3";

	readAllFingerprints(database);
	makeAllVaults("triangles", 1.0);

	GenProps gp = nextGeneration(newGeneration(populationSize, genomeSize), populationSize, genomeSize, numFingerprints, distroIncrement, database);

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
	cout<<"readAllFingerprintsTest:\t\t\t";
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
				readAllFingerprints(goodDatabase);
				vector<vector<vector<Minutiae> > > goodFingerprints = allFingerprints;
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
		readAllFingerprints(badDatabase);
		vector<vector<vector<Minutiae> > > badFingerprints = allFingerprints;
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
	cout<<"allFingerprintsReadCorrectlyTest:\t\t";
	cout.flush();
	bool toReturn = true;

	readAllFingerprints("FVC2002/DB1");
	makeAllVaults("triangles", 1.0);

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

	readAllFingerprints("FVC2002/DB1");
	makeAllVaults("triangles", 1.0);

	toReturn = toReturn && (allVaults.size()==100);
	for(int i = 0; i<allVaults.size(); i++)
	{
		toReturn = toReturn && (allVaults.at(i).size()==8);
	}
	for(int i = 0; i<allVaults.size(); i++)
	{
		for(int j = 0; j<allVaults.at(i).size(); j++)
		{
			toReturn = toReturn && ((allVaults.at(i).at(j)->unlock(allFingerprints.at(i).at(j))).score>0);
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
	cout<<"runDatabaseTest:\t\t\t\t";
	cout.flush();
	bool toReturn = true;

	int genomeSize = 20;
	string database = "FVC2004/DB4";
	int numFingerprints = 50;
	float distroIncrement = 1.1;
	vector<float> distro = newDistro(genomeSize);

	readAllFingerprints("FVC2002/DB1");
	makeAllVaults("triangles", 1.0);

	float eer = runDatabase(database, &distro, numFingerprints, distroIncrement);
	toReturn = toReturn && (eer < 1) && (eer>0);
	if(toReturn)
		cout<<"passed"<<endl;
	else
		cout<<"failed"<<endl;

	return toReturn;
}


void runGAUnitTests(int *count, int *pass){

	(*count)++; *pass += evaluatePerformanceTest();
	(*count)++; *pass += nextGenerationTest();
	(*count)++; *pass += readAllFingerprintsTest();
	(*count)++; *pass += allFingerprintsReadCorrectlyTest();
	(*count)++; *pass += allVaultsReadCorrectlyTest();
	(*count)++; *pass += runDatabaseTest();


}
