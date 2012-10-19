//./RUNME genetic 1 20 1 1 FVC2002/DB1 10 bf

#include <iostream>
#include <stdlib.h>
#include "../common/inc/RunTest.h"
#include "../common/inc/ProbDistroGA.h"
#include "../common/inc/UnitTests.h"
#include "../common/inc/IOUtils.h"
#include "../matching_tjea/inc/matching.h"
#include "../matching_tjea/inc/global_parameters.h"

using namespace std;

/*
currently prints to stdout, and you must find the best
threshold by inspection
*/
void findBestThreshold(string methodType, int genSize=200, float distroIncrement=.25, string database="FVC2002/DB1")
{
	//@TODO file io
	//initialize the NTL stuff
	ZZ_p::init(FEILD_SIZE);
	srand(time(NULL));
	SetSeed(to_ZZ(rand()));

	int popSize = 1, iterations = 1, numFingerprints=100;
	vector<float> thresholds = vector<float>(genSize);
	makeRoot();
	readAllFingerprints(database);
	makeAllFolded(methodType);
	//cout << endl << "start." << endl << flush;
	makeAllFoldedTransformed(methodType);
	//cout << "END." << endl << flush;
	makeAllVaults(methodType, distroIncrement);
	float bestSoFar = 1.0;
	vector<float> wtf = vector<float>(10000);

	for(int i = 0; i<genSize; i++)
	{
		thresholds.at(i) = 1;
		cout<<"thresh = "<<i/(1/distroIncrement)<<" eer =\t"<<flush;
		float currentEER = 1/evalFitness(thresholds, numFingerprints, distroIncrement, database);
		cout<<currentEER;
		if(currentEER<bestSoFar)
		{
			bestSoFar=currentEER;
			cout<<"\t<< best so far!";
		}
		cout<<endl;
	}
	cout<<endl<<"done."<<endl;
}

void printHelp()
{
	cout<<endl;
	cout<<"to print these instructions, type:"<<endl;
	cout<<"\t[exec] \"help\" "<<endl<<endl;
	cout<<"to run a genetic algorithm, type:"<<endl;
	cout<<"\t[exec] \"genetic\", int populationSize, int genomeSize, float distroIncrement, int iterations, "
				<<endl<<"\t\t\tstring pathToDB, int numFingerprints, string method"<<endl;
	cout<<"\te.g. ./RUNME genetic 2 20 1.0 10 FVC2002/DB1 10 triangles"<<endl;
	cout<<"to use the default parameters, type"<<endl
				<<"\t[exec] genetic default"<<endl<<endl;
	cout<<"to run an experiment to find the best threshold, type:"<<endl;
	cout<<"\t[exec] \"findBestThreshold\" int genSize, float distroIncrement, string database, string methodType"<<endl;
	cout<<"\te.g. ./RUNME findBestThreshold 10 .5 FVC2002/DB1 bf"<<endl;
	cout<<"to use the default parameters, type"<<endl
				<<"\te.g. [exec] findBestThreshold default methodType"<<endl<<endl;

	cout<<"method type is one of {bf, triangles}"<<endl<<endl;
}

void runGenetic(int populationSize, int genomeSize, float distroIncrement, int iterations, 
				string database, int numFingerprints, string methodType, string call)
{
	makeRoot();
	readAllFingerprints(database);
	//makeAllRandomFingerprints();
	//cout << endl << "start." << endl << flush;
	makeAllFolded(methodType);
	//makeAllFoldedTransformed(methodType);
	makeAllFoldedRotated(methodType);
	//cout << "END." << endl << flush;
	makeAllVaults(methodType, distroIncrement);
	cout
			//<<endl
			<<"populationSize = "<<populationSize<<endl
			<<"genomeSize = "<<genomeSize<<endl
			<<"distroIncrement = "<<distroIncrement<<endl
			<<"iterations = "<<iterations<<endl
			<<"database = "<<database<<endl
			<<"numFingerprints = "<<numFingerprints<<endl
			<<"method = "<<methodType<<endl<<endl;
/*
struct FPTemplate
{
	Minutiae marr[CUBS_MAX_MINUTIAE];
	int		 nCnt;
};
typedef struct FPTemplate* ptrFPTemplate;


int CUBS_MatchFeaturesInternal(ptrFPTemplate hRef, ptrFPTemplate hTst, 
						         	   MatchResult*	pMatchResult, MatchResultEx* pResultEx);
*/



	beginGenetic(populationSize, genomeSize, distroIncrement, iterations, database, numFingerprints, call);
}


//type [exec] "help" for more information
int main(int argc, char** argv)
{

	//initialize the NTL stuff
	ZZ_p::init(FEILD_SIZE);
	srand(time(NULL));
	SetSeed(to_ZZ(rand()));

	Global_Parameters G_P();

	if(argc==2 && string(argv[1]) == "help")
	{
		printHelp();
		return 0;
	}
	else if(argc==2 && string(argv[1]) == "test")
	{
		runUnitTests();
		return 0;
	}
	else if(argc >= 4 && string(argv[1]) == "findBestThreshold")
	{
		if(argc == 4 && string(argv[2]) == "default")
			findBestThreshold(string(argv[3]));
		else if (argc == 6)
		{
			findBestThreshold(argv[5], atoi(argv[2]), atof(argv[3]), argv[4]);
		}
		else
		{
			cout	<<"command not understood. type"<<endl<<"\t[exec] help "<<endl
					<<"to see the command line options to this program"<<endl;
		}
		return 0;
	}
	else if(argc >= 2 && string(argv[1]) == "genetic")
	{
		int populationSize, genomeSize, iterations, numFingerprints;
		float distroIncrement;
		string database, methodType;

		if(argc == 9)
		{
			populationSize = atoi(argv[2]);
			genomeSize = atoi(argv[3]);
			distroIncrement = atof(argv[4]);
			iterations = atoi(argv[5]);
			database = argv[6];
			numFingerprints = atoi(argv[7]);
			methodType = argv[8];
			runGenetic(populationSize, genomeSize, distroIncrement, iterations, database, numFingerprints, methodType, concatenateNotFirst(argc, argv));
		}

		else if (argc == 3 && string(argv[2])=="default")
		{
			cout<<"using defaults"<<endl;
			runGenetic(2, 20, 1.0, 10, "FVC2002/DB1", 10, "triangles", "genetic_default");
		}
		else
		{
			cout	<<"command not understood. type"<<endl<<"\t[exec] help "<<endl
					<<"to see the command line options to this program"<<endl;
		}
		return 0;
	}

	else
	{
		cout<<"type ./RUNME help for more instructions on experiments you can run"<<endl<<endl;
		return 0;
	}

	cout<<"done."<<endl;
}

