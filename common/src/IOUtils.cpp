#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>
//so that we can read the database of fingerprints into a static variable
#include "../inc/RunTest.h"
#include "../inc/Common.h"
#include "../inc/fingerprint_model.h"
#include "../inc/Vault.h"

using namespace std;

/*used to remember what the program's args were
so it omits the first entry (argv is passed as arr)
*/
string concatenateNotFirst(int len, char* arr[])
{
   if (len < 1) {
       return "";
   }
   string result = "";
   for (int i=1; i < len-1; ++i) {
       result += arr[i];
       result += "_";
   }
	result+=arr[len-1];

	for(int i = 1; i<result.length(); i++)
		if (result.at(i)=='/')
			result.at(i) = '_';
	
    return result;
}

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

void print(minutia m)
{
	cout<<m.x<<", "<<m.y<<", "<<m.theta<<", "<<endl;
}

void print(vector<minutia> v)
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

vector<minutia> readFingerprintFromFile(string inputFilename)
{
	try
	{
		minutia	currentMinutiae;
		vector<minutia> toReturn;
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
					currentMinutiae.x = atoi(tokenized.at(0).c_str());
					currentMinutiae.y = atoi(tokenized.at(1).c_str());
					currentMinutiae.theta = atoi(tokenized.at(2).c_str());
				}

				toReturn.push_back(currentMinutiae);
			}
			inputFile.close();
		}

		//This is because the first 5 lines that get read from each file are trash
		return vector<minutia>(toReturn.begin()+5, toReturn.end());
	}
	catch(bad_alloc& ba)
	{
		//cerr<<inputFilename<<endl;
		return vector<minutia>();
	}
}

//@TODO references
//@param database is like "FVC2002/DB1"
void readAllFingerprints(string database/*, vector<vector<vector<minutia> > > allFingerprints*/)
{
	vector<vector<vector<minutia> > > allFP;
	for(int fingerprintNumber=1; fingerprintNumber<=100; fingerprintNumber++)
	{
		vector<vector<minutia> > toPush;
		for(int readingNumber=1; readingNumber<=8; readingNumber++)
		{
			vector<minutia> reddFingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/" 
				+ database +"/features/" + toString(fingerprintNumber) + "_" + 
				toString(readingNumber) + ".fp");
			//cout << endl << "before: " << reddFingerprint[0].x << " " << reddFingerprint[0].y;
			center(reddFingerprint);
			//cout << endl << "after: " << reddFingerprint[0].x << " " << reddFingerprint[0].y;
			toPush.push_back(reddFingerprint);
		}
		allFP.push_back(toPush);
	}
	allFingerprints = allFP;
}


void center(vector<minutia>& print)
{
	int center_x = FP_IMAGE_WIDTH>>1;
	int center_y = FP_IMAGE_HEIGHT>>1;
	int n=print.size();
	int sumx = 0;
	int sumy = 0;

	for(int i=0; i<n; i++){
		sumx += print[i].x;
		sumy += print[i].y;
	}

	int averagex = sumx/n; 
	int averagey = sumy/n;
	//cout << averagex << " " << averagey << endl;

	int diffx = averagex - center_x;
	int diffy = averagey - center_y;
	
	for(int i=0; i<n; i++){
		print[i].x -= diffx;
		print[i].y -= diffy;
	}

} 


void makeAllRandomFingerprints(/*string database, vector<vector<vector<minutia> > > allFingerprints*/)
{
	vector<vector<vector<minutia> > > allFP;
	FingerprintModel model;
	for(int fingerprintNumber=1; fingerprintNumber<=100; fingerprintNumber++)
	{
		vector<vector<minutia> > toPush;
		model.GenerateRandom();
		for(int readingNumber=1; readingNumber<=8; readingNumber++)
		{
			vector<minutia> reddFingerprint = vector<minutia>();
			model.MakePrint(reddFingerprint);
			toPush.push_back(reddFingerprint);
			//print(reddFingerprint);
			
			//cout << reddFingerprint.size() << endl;
		}
		allFP.push_back(toPush);
	}
	allFingerprints = allFP;
}
