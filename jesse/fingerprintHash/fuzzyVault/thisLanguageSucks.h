#ifndef __THISLANGUAGESUCKS_H__
#define __THISLANGUAGESUCKS_H__


#include <vector>
#include <algorithm>
#include <fstream>
#include <set>
#include <time.h>
#include <limits>

#include "NTL/mat_ZZ.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

typedef struct
{
	ZZ_p z1;
//	ZZ_p z2;
	ZZ_p f1;
//	ZZ_p f2;
	bool chaff;
}folded;

typedef struct
{
	vec_ZZ_p poly;
	double score;
}thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture;
/*
typedef struct
{
	int		m_nFeatureNo;   //index of this feature
    int     m_nX;           //x -coordinate
    int     m_nY;           //y -coordinate
    int     m_nTheta;       //orientation
    char    m_cType;        //minutiae type
    int     m_nQuality;     //minutiae quality
}Minutiae;  // Used in CUBSEnroll and CUBSMatch

const int CUBS_MAX_MINUTIAE = 255;

struct FPTemplate
{
	Minutiae marr[CUBS_MAX_MINUTIAE];
	int		 nCnt;
};

typedef struct FPTemplate* ptrFPTemplate;



typedef struct
{
  bool   MatchPerformed;
  double similarity;
}MatchResult;

typedef struct
{
	int		Lidx; // feature number in first list
	int		Ridx; // feature number in second list
	double	score; // matching score
}FeaPair;

typedef struct
{
	int     NumOfMatchedMinutiae;
	double	DistScore;
	double	PointScore;
	int		NumOfOverlap_a; // number of overlapped minutiae on A
	int		NumOfOverlap_b; // number of overlapped minutiae on B
	int		NumOfMinutiae_a;  // number of minutiae on A
	int		NumOfMinutiae_b;  // number of minutiae on B
	int     Width_a;          // width of print A (input)
	int     Height_a;         // height of print A (input)
	int     Width_b;          // width of print B (reference)
	int     Height_b;         // height of print B (reference)
	int     Width_comb;       // width of combined print
	int     Height_comb;      // height of combined print
	FeaPair FeaCorr[CUBS_MAX_MINUTIAE]; // Feature correspondence 
}MatchResultEx;
*/

#include "common.h"
#include "WelchBerlekamp.h"
#include "math.h"
//#include "global_parameters.h"
#include "VaultMethod.h"


using namespace std;






#endif