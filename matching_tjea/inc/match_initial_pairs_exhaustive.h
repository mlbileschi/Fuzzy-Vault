#ifndef __MATCH_INITIAL_PAIRS_H__
#define __MATCH_INITIAL_PAIRS_H__

#include "../inc/match_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

#define MAX_R 100	 /*points farther away from this are not matched*/

/*!
	@function UTLConvertToRadial
	@abstract Converts the minutiae array from cartesian to radial
*/
int	UTLConvertToRadial(const Minutia pMin[], const int nCnt,const int nX,const int nY, const int nTheta, RMinutia pRMin[]);

/*!
	@function UTLMatchRadial
	@abstract Exhaustively matches all possible pairs of radial minutia and popluates the pair matrix
*/
int UTLMatchRadial(const RMinutia pRMin[], const int nRCnt, const RMinutia pTMin[], const int nTCnt,Pair pPair[]);

/*!
	@function UTLSortPairs
	@abstract sorts the pair scores
*/
int UTLSortPairs(Pair pPair[], const int nCnt);


#if defined(__cplusplus)
}
#endif

#endif
