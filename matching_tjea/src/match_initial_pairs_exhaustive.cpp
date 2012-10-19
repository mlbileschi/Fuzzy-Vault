#include <stdlib.h>
#include <math.h>

#include "../inc/match_types.h"
#include "../inc/match_utils.h"
#include "../inc/match_initial_pairs_exhaustive.h"

/*-------------------------------------------------------------------------------------------------------------
UTLConvertToRadial
-------------------------------------------------------------------------------------------------------------*/
int	UTLConvertToRadial(const Minutia pMin[], const int nCnt,const int nX,const int nY, const int nTheta, RMinutia pRMin[])
{
	int				i	=	0;
	const double	PI	=	3.14159;
	/*sanity check*/
	if(pMin == NULL || nCnt < 1 || nX < 0 || nY < 0 || nTheta < 0 || pRMin == NULL)
		return UTL_ERR_INVALID_PARAMETER;

	for(i=0;i<nCnt;i++)
	{
		pRMin[i].m_dR		=	sqrt((pMin[i].m_nX-nX)*(pMin[i].m_nX-nX)+(pMin[i].m_nY-nY)*(pMin[i].m_nY-nY));
		pRMin[i].m_nTheta	=	(int)(atan2(pMin[i].m_nY-nY,pMin[i].m_nX-nX)*180/PI);
		if(pRMin[i].m_nTheta < 0)
		{
			pRMin[i].m_nTheta += 360;
		}
		pRMin[i].m_nDelta	=	pMin[i].m_nTheta-nTheta;
		if(pRMin[i].m_nDelta < 0)
		{
			pRMin[i].m_nDelta += 360;
		}
	}
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLConvertToRadial
-------------------------------------------------------------------------------------------------------------*/
int UTLMatchRadial(const RMinutia pRMin[], const int nRCnt, const RMinutia pTMin[], const int nTCnt,Pair pPair[])
{
	int i			=	0;
	int j			=	0;
	double  dR		=	0;
	double	dTheta	=	0;
	double	dDelta	=	0;
	int		nAbs	=	0;
	int		nTFlags[MAX_MINUTIAE]	={0};
	int		nRFlags[MAX_MINUTIAE]	={0};


	for(i=0;i<nRCnt;i++)
	{
		/*skip far away points*/
		if(pRMin[i].m_dR > MAX_R)
			continue;

		for(j=0;j<nTCnt;j++)
		{
			if(nTFlags[j]) /*already matched*/
				continue; 

			/*skip far away points*/
			if(pTMin[j].m_dR > MAX_R)
				continue;

			/*check for match*/
			dR		=	fabs(log10((pRMin[i].m_dR+1e-5)/(pTMin[j].m_dR+1e-5)));
				
			if(dR > DR) /*no match*/
				continue;

			nAbs	=	abs(pRMin[i].m_nTheta - pTMin[j].m_nTheta);
			if((360-nAbs) < nAbs)
				dTheta = (360-nAbs);
			else
				dTheta = nAbs;

			if(dTheta > DTHETA) /*no match*/
				continue;

			nAbs	=	abs(pRMin[i].m_nDelta - pTMin[j].m_nDelta);
			if((360-nAbs) < nAbs)
				dDelta = (360-nAbs);
			else
				dDelta = nAbs;

			if(dDelta > DDELTA) /*no match*/
				continue;
			
			
			if(!nRFlags[i] && !nTFlags[j]) /*match found*/
			{
				nRFlags[i]	=	1;
				nTFlags[j]	=	1;
				//add vote to pair matrix
				pPair[i*nTCnt+j].m_dScore++;
				pPair[i*nTCnt+j].m_nRIdx	=	i;
				pPair[i*nTCnt+j].m_nTIdx	=	j;
				DBGPrintf("Radial Match: [%d,%d]->(%lf,%d,%d),(%lf,%d,%d)\n",i,j,pRMin[i].m_dR,pRMin[i].m_nTheta,pRMin[i].m_nDelta,
							pTMin[j].m_dR,pTMin[j].m_nTheta,pTMin[j].m_nDelta);
			}
		}
	}
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLPairSortingFunction
-------------------------------------------------------------------------------------------------------------*/
int UTLPairSortingFunction(const void* p,const void* q)
{
	const Pair* pPairA = (const Pair*)p;
	const Pair* pPairB = (const Pair*)q;
	if(pPairA->m_dScore  < pPairB->m_dScore)
		return 1;
	else
		return -1; //descending order
}

/*-------------------------------------------------------------------------------------------------------------
UTLSortPairs
-------------------------------------------------------------------------------------------------------------*/
int UTLSortPairs(Pair pPair[], int nCnt)
{
	/*sanity check*/
	if(pPair == NULL || nCnt < 1) 
		return UTL_ERR_INVALID_PARAMETER;
	qsort(pPair,nCnt,sizeof(Pair),UTLPairSortingFunction);
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLGetInitialPairs
-------------------------------------------------------------------------------------------------------------*/
int UTLGetInitialPairs(const Minutia pminR[], const int nRCnt, const Minutia pminT[], const int nTCnt, Pair pPairs[])
{
	int		i	=	0;
	int		j	=	0;
	
	RMinutia rminR[MAX_MINUTIAE];
	RMinutia rminT[MAX_MINUTIAE];

	/*sanity check*/
	if(pminR == NULL || nRCnt < 1 || pminT == NULL || nTCnt < 1 || pPairs == NULL)
		return UTL_ERR_INVALID_PARAMETER;
	
	/* test matching hypotheses */
	for(i=0;i<nRCnt;i++)
	{
		UTLConvertToRadial(pminR,nRCnt,pminR[i].m_nX,pminR[i].m_nY,pminR[i].m_nTheta,rminR);
		for(j=0;j<nTCnt;j++)
		{
			UTLConvertToRadial(pminT,nTCnt,pminT[j].m_nX,pminT[j].m_nY,pminT[j].m_nTheta,rminT);
			DBGPrintf("----------------------------------------------------\n");
			DBGPrintf("MatchMinutiaSets: Hypothesis (%d,%d)->(%d,%d)\n",pminR[i].m_nX,pminR[i].m_nY,pminT[j].m_nX,pminT[j].m_nY);
			DBGPrintf("----------------------------------------------------\n");
			UTLMatchRadial(rminR,nRCnt,rminT,nRCnt,pPairs);
		}
	}
	UTLSortPairs(pPairs,nRCnt*nTCnt);
	DBGPrintPairs(pminR,nRCnt,pminT,nRCnt,pPairs);
	return UTL_ERR_NO_ERROR;
}
