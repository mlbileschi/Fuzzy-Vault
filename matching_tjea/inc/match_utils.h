#ifndef __MATCH_UTILS_H__
#define __MATCH_UTILS_H__

#include "../inc/match_types.h"




#if defined(__cplusplus)
extern "C" {
#endif


/*!
	@function  UTLCreateMinutiaList
	@abstract  Creates a list of minutiae from array of co-ordinates and angles
*/
int	UTLCreateMinutiaList(const int pnX[],const int pnY[], const int pnTheta[], const int nCnt, Minutia pMin[]);

/*!
	@function UTLGetTransformationParams
	@abstract Gets the transformation parameters
*/
int UTLGetTransformationParams(const Minutia pMinR[], const Minutia pMinT[],const Pair pPair[], const int nCnt, TxParams* pTx);

/*!
	@function UTLTransformMinutiaSet
	@abstract transforms the minutia based on the TPS parameters
*/
int UTLTransformMinutiaSet(Minutia pMin[],const int nCnt,const TxParams* pTx);

/*!
	@function UTLFindMatches
	@abstract updates the pair matrix with new matches
*/
int UTLFindMatches(Minutia pMinR[],int nRCnt,Minutia pMinT[],int nTCnt,Pair pPairs[],int* pnCnt);


/*!
	@function UTLGetMatchScore
	@abstract Compares the two set of minutiae for matc
*/
int UTLGetMatchScore(const Minutia pRMin[],const int nRCnt,const Minutia pTMin[], const int nTCnt,double* pdScore);

/*!
	@function UTLGetMatchScoreRadial
	@abstract Compares the two set of radial minutia
*/
int UTLGetMatchScoreRadial(const RMinutia pRMin[],const int nRCnt,const RMinutia pTMin[],const int nTCnt,double* pdScore);

/*!
	@function UTLGetInitialPairs
	@abstract Gets the potential corresponding pairs
*/
int UTLGetInitialPairs(const Minutia pMinR[], const int nRCnt, const Minutia pMinT[], const int nTCnt, Pair pPairs[]);


/*!
	@function MatchMinutiaSets
	@abstract Driver routine
*/
int MatchMinutiaSets(const int pnRX[], const int pnRY[], const int pnRTheta[], int nRCnt,
					 const int pnTX[], const int pnTY[], const int pnTTheta[], int nTCnt,
					 double *pdScore);


/*!
	@function UTLSortPairs
	@abstract 
*/
int UTLSortPairs(Pair pPair[], int nCnt);

#if defined(__cplusplus)
}
#endif

#endif
