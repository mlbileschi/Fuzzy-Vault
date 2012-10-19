#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>
//so that we can read the database of fingerprints into a static variable
#include "../inc/runTest.h"
#include "../inc/common.h"

using namespace std;

void print(vector<float> v)
{
	for( int i = 0; i<v.size(); i++)
	{
		cout<<v.at(i)<<", ";
	}
	cout<<endl<<endl;
}

void print(vector<string> v)
{
	for( int i = 0; i<v.size(); i++)
	{
		cout<<v.at(i)<<", ";
	}
	cout<<endl;
}

void print(Minutiae m)
{
	cout<<m.m_nX<<", "<<m.m_nY<<", "<<m.m_nTheta<<", "<<endl;
}

void print(vector<Minutiae> v)
{
	for(int i = 0; i<v.size(); i++)
	{
		print(v.at(i));
		cout<<endl;
	}
}

string toString(int number)
{
   stringstream ss;
   ss << number; 
   return ss.str();
}

vector<string> splitBySpaces(string toSplit)
{
	string toPush;
	vector<string> toReturn;
	istringstream toSplitISS(toSplit);
	while(getline(toSplitISS, toPush, ' '))
	{
		toReturn.push_back(toPush);
	}
	return toReturn;
}

vector<Minutiae> readFingerprintFromFile(string inputFilename)
{
	try
	{
		Minutiae	currentMinutiae;
		vector<Minutiae> toReturn;
		string line;
		ifstream inputFile (inputFilename.c_str());
		if (inputFile.is_open())
		{
			while ( inputFile.good() )
			{
				int position = 0;
				getline (inputFile, line);
				vector<string> tokenized = splitBySpaces(line);
				if(tokenized.size()==5)
				{
					currentMinutiae.m_nX = atoi(tokenized.at(0).c_str());
					currentMinutiae.m_nY = atoi(tokenized.at(1).c_str());
					currentMinutiae.m_nTheta = atoi(tokenized.at(2).c_str());
				}

				toReturn.push_back(currentMinutiae);
			}
			inputFile.close();
		}

		//This is because the first 5 lines that get read from each file are trash
		return vector<Minutiae>(toReturn.begin()+5, toReturn.end());
	}
	catch(bad_alloc& ba)
	{
		//cerr<<inputFilename<<endl;
		return vector<Minutiae>();
	}
}

//@TODO references
//@param database is like "FVC2002/DB1"
void readAllFingerprints(string database/*, vector<vector<vector<Minutiae> > > allFingerprints*/)
{
	vector<vector<vector<Minutiae> > > allFP;
	for(int fingerprintNumber=1; fingerprintNumber<=100; fingerprintNumber++)
	{
		vector<vector<Minutiae> > toPush;
		for(int readingNumber=1; readingNumber<=8; readingNumber++)
		{
			vector<Minutiae> reddFingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/" +
																		database +"/features/" + 
																		toString(fingerprintNumber) + "_" + 
																		toString(readingNumber) + ".fp");
			toPush.push_back(reddFingerprint);
		}
		allFP.push_back(toPush);
	}
	allFingerprints = allFP;
}
