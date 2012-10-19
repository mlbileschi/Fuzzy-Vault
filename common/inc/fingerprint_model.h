#ifndef FINGERPRINT_MODEL_H
#define FINGERPRINT_MODEL_H

#include <vector>
//#include "common.h"
#include "../inc/Common.h"



class FingerprintModel
{
public:
	int GenerateRandom();
	int MakePrint(vector<minutia>& fp);

	vector<minutia> original_finger;
};



#endif // FINGERPRINT_MODEL_H

