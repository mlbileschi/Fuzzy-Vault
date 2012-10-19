#include <vector>
#include <algorithm>
#include "../inc/evaluatePerformance.h"
#include "../inc/IOUtils.h"

using namespace std;

int is_this_working(){return 15;}

float evaluatePerformance(vector<float> impostors, vector<float> genuines)
{

	sort (genuines.begin(), genuines.end());
	sort (impostors.begin(), impostors.end());

	float genInc = 1.0/genuines.size(); //genuine increment. applies to FRR
	float impDec = 1.0/impostors.size(); //impostor decrement. applies to FAR
	int genPtr = 0;
	int impPtr = 0;
	float curFAR = 1.0;
	float curFRR = 0.0;
	vector<float> FAR;
	vector<float> FRR;

	while(true)
	{
		if(impPtr==impostors.size()-1 && genPtr==genuines.size()-1) //we have reached the end
		{
			break;
		}
		else if(impPtr==impostors.size()-1 && genPtr!=genuines.size()-1) //work on genuines
		{
			curFRR+=genInc;
			genPtr++;
			FAR.push_back(curFAR);
			FRR.push_back(curFRR);
			continue;			
		}
		else if(impPtr!=impostors.size()-1 && genPtr==genuines.size()-1) //work on genuines
		{
			curFAR-=impDec;
			impPtr++;
			FRR.push_back(curFRR);
			FAR.push_back(curFAR);
			continue;
		}
		else //we may increment either, and we should check which one we should do
		{
			if(impostors.at(impPtr) < genuines.at(genPtr))
			{
				curFAR-=impDec;
				impPtr++;
				FAR.push_back(curFAR);
				FRR.push_back(curFRR);
				continue;
			}
			else if(impostors.at(impPtr) > genuines.at(genPtr))
			{
				curFRR+=genInc;
				genPtr++;
				FAR.push_back(curFAR);
				FRR.push_back(curFRR);
				continue;
			}
			else //they are equal so do both
			{
				curFAR-=impDec;
				impPtr++;
				FAR.push_back(curFAR);
				FRR.push_back(curFRR);
				curFRR+=genInc;
				genPtr++;
				FAR.push_back(curFAR);//
				FRR.push_back(curFRR);
				continue;				
			}
		}		
	}

	float EER = 999;
	for(int i = 0; i < min(FAR.size(), FRR.size()); i++)
	{

		//cout << FAR.at(i) << " \t" << FRR.at(i) << endl;

		if(FAR.at(i)<FRR.at(i))
		{
			EER = (FAR.at(i) + FRR.at(i))/2.0;
			break;
		}
		else if (FAR.at(i)==FRR.at(i))
		{
			EER = FAR.at(i);
			break;
		}
	}

	return EER;
}
