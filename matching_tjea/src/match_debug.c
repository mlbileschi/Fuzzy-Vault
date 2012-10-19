#include "../inc/match_types.h"
#include "../inc/match_debug.h"
#include <stdarg.h>
#include <stdio.h>

#if MATCH_DEBUG

//--------------------------------------------------------------
//DBGPrintf
//--------------------------------------------------------------
void DBGPrintf(const char* pcFormat,...) /*will be different for unix */
{
	FILE* fp = fopen("debug.txt","wa");
	va_list		argptr;
	va_start(argptr,pcFormat);
	//vprintf(pcFormat,argptr);
	vfprintf(fp,pcFormat,argptr);
	fclose(fp);
	va_end(argptr);
}
//--------------------------------------------------------------
//DBGPrintPairs
//--------------------------------------------------------------
void DBGPrintPairs(Minutia pMinR[],int nRCnt,Minutia pMinT[],int nTCnt,Pair pairs[])
{
	int i;
	DBGPrintf("---------------------------------------------\n");
	DBGPrintf("Pair info\n");
	DBGPrintf("---------------------------------------------\n");
	for(i=0;i<nRCnt*nTCnt;i++)
	{
		DBGPrintf("[%d,%d]->(%d,%d)-(%d,%d)-->%lf\n",pairs[i].m_nRIdx,pairs[i].m_nTIdx,pMinR[pairs[i].m_nRIdx].m_nX,
											pMinR[pairs[i].m_nRIdx].m_nY,pMinT[pairs[i].m_nTIdx].m_nX,
											pMinT[pairs[i].m_nTIdx].m_nY,pairs[i].m_dScore);
	}
}
//--------------------------------------------------------------
//DBGPrintMinutia
//--------------------------------------------------------------
void DBGPrintMinutia(const Minutia pMin[],const int nCnt)
{
	int i=0;
	DBGPrintf("-----------------------------------------------\n");
	DBGPrintf("Minutia details\n");
	DBGPrintf("-----------------------------------------------\n");
	for(i=0;i<nCnt;i++)
	{
		DBGPrintf("(%d,%d,%d)\n",pMin[i].m_nX,pMin[i].m_nY,pMin[i].m_nTheta);
	}
}

#else
//--------------------------------------------------------------
//DBGPrintf
//--------------------------------------------------------------
void DBGPrintf(const char* pFormat,...) /*will be different for unix */
{}

//--------------------------------------------------------------
//DBGPrintPairs
//--------------------------------------------------------------
void DBGPrintPairs(Minutia pMinR[],int nRCnt,Minutia pMinT[],int nTCnt,Pair pairs[])
{}

//--------------------------------------------------------------
//DBGPrintMinutia
//--------------------------------------------------------------
void DBGPrintMinutia(const Minutia pMin[],const int nCnt)
{}

#endif
