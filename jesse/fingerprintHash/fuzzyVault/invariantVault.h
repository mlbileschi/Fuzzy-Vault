#pragma once

#ifndef __INVARIANTVAULT_H__
#define __INVARIANTVAULT_H__

#include "vault.h"

/*
typedef struct
{
	double d;
	double sigma;
}dsigma;
*/

typedef struct
{
	double d[K];
	double sigma[K];
//	dsigma dsigma[K];
	ZZ_p f1;
	ZZ_p f2;
	bool chaff;
}triangles;


void lockInv(vector<triangles>* vault, ptrFPTemplate hRef, vec_ZZ_p secret);
thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture unlockInv(vector<triangles>* vault, ptrFPTemplate hTst);

bool comparer(vector<triangles> f1, vector<triangles> f2);

vector<triangles> Minutiae2triangles(ZZ x, ZZ y, ZZ theta, ptrFPTemplate temp);
ZZ_p Triangles2zzp(triangles tri);

#endif //__INVARIANTVAULT_H__