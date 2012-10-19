#include <iostream>
#include "../common/inc/runTest.h"
#include "../common/inc/ProbDistroGA.h"
#include "../common/inc/unitTests.h"


/* This was run with a broken setDistro so the best distro is irrelevant. Used triangles

Jesses-MacBook-Pro:main jessehartloff$ ./RUNME 5 20 .5 200 FVC2002/DB2 20

populationSize = 5
genomeSize = 20
distroIncrement = 0.5
iterations = 200
database = FVC2002/DB2
numFingerprints = 20

9 7 5 17 15 BEGIN GENETIC ALGORITHM
fitness of best = 0.188417	convergence : 172.314	  99% finished
END GENETIC ALGORITHM

best was: 
(0.916422,0.811337,0.779064,0.763556,0.752269,0.728852,0.622558,0.615614,0.605236,0.553948,0.451753,0.393037,0.33437,0.311365,0.202809,0.198493,0.142793,0.109702,0.0914791,0.0663683,)
with score : 0.188417

done.
*/

using namespace std;

//used to remember what the program's args were
string concatenate(int argc, char* argv[])
{
   if (argc < 1) {
       return "";
   }
   string result = "";
   for (int i=1; i < argc-1; ++i) {
       result += argv[i];
       result += "_";
   }
	result+=argv[argc-1];

	for(int i = 1; i<result.length(); i++)
		if (result.at(i)=='/')
			result.at(i) = '_';
	
    return result;
}



//sysargs: int populationSize, int genomeSize, float distroIncrement, int iterations, 
//				string database, int numFingerprints
int main(int argc, char* argv[])
{

	//initialize the NTL stuff
	ZZ_p::init(FEILD_SIZE);
	srand(time(NULL));
	SetSeed(to_ZZ(rand()));

	if(argc==2 && string(argv[1]) == "test")
	{
		runUnitTests();
		return 1;
	}


	if(argc!=8)
	{
		cout<<endl<<"usage: [exec] int populationSize, int genomeSize, float distroIncrement, int iterations, "
					<<"string database, int numFingerprints, method"<<endl<<endl;
		cout<<"e.g. ./RUNME 2 20 1.0 10 FVC2002/DB1 10 triangles"<<endl<<endl;
		return 1;
	}

	readAllFingerprints(argv[5]);
	makeAllFolded(argv[7]);
cout << endl << "start." << endl << flush;
	makeAllFoldedTransformed(argv[7]);
int crap;
cout << "END." << "    size of folded, z1: " << sizeof(allFingerprintsFolded[0][0][0]) << "," << sizeof(allFingerprintsFolded[0][0][0].z1) << sizeof(allFingerprintsFolded[0][0][0].f1) << sizeof(allFingerprintsFolded[0][0][0].chaff)  << ", " << sizeof(crap) << endl << flush;
	makeAllVaults(argv[7], atof(argv[3]));
	cout<<endl
			<<"populationSize = "<<atoi(argv[1])<<endl
			<<"genomeSize = "<<atoi(argv[2])<<endl
			<<"distroIncrement = "<<atof(argv[3])<<endl
			<<"iterations = "<<atoi(argv[4])<<endl
			<<"database = "<<argv[5]<<endl
			<<"numFingerprints = "<<atoi(argv[6])<<endl
			<<"method = "<<argv[7]<<endl<<endl;

	begin(atoi(argv[1]), atoi(argv[2]), atof(argv[3]), atoi(argv[4]), argv[5], atoi(argv[6]), concatenate(argc, argv));

	cout<<endl<<"done."<<endl;
	return 0;
}
