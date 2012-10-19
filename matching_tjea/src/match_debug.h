#ifndef __MATCH_TYPES_H__
#define __MATCH_TYPES_H__


#if defined(__cplusplus)
extern "C" {
#endif

/*!
	@function DBGPrintf
	@abstract debug print
*/
void DBGPrintf(const char*,...);
void DBGPrintPairs(Minutia pMinR[],int nRCnt,Minutia pMinT[],int nTCnt,Pairs pairs[]);
void DBGPrintMinutia(const Minutia pMin[],const int nCnt);
#if defined(__cplusplus)
}
#endif


#endif