#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>

#include "../inc/match_types.h"
#include "../inc/match_utils.h"
#include "../inc/match_debug.h"
#include "../inc/f2c.h"
#include "../inc/lapack.h"

/*-------------------------------------------------------------------------------------------------------------
UTLCreateMinutiaList
-------------------------------------------------------------------------------------------------------------*/
int	UTLCreateMinutiaList(const int pnX[],
						 const int pnY[], 
						 const int pnTheta[], 
						 const int nCnt,
						 Minutia pMin[])
{
	int		i	=	0;
	/*sanity check*/
	if(pnY == NULL || pnY == NULL || pnTheta == NULL || nCnt < 1 || pMin == NULL)
		return UTL_ERR_INVALID_PARAMETER;

	for(i=0;i<nCnt;i++)
	{
		pMin[i].m_nX	=	pnX[i];
		pMin[i].m_nY	=	pnY[i];
		pMin[i].m_nTheta=	pnTheta[i];
	}
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLGetTransformationParams
-------------------------------------------------------------------------------------------------------------*/
int UTLGetTransformationParams(const Minutia pMinR[], const Minutia pMinT[],const Pair pPair[], const int nPairCnt, TxParams* pTx)
{
	int			nRow							=	0;
	int			nCol							=	0;
	const int	nMax							=	MAX_TPS_ORDER-3;
				
	int			nOrder							=	0;
	int			idx								=	0;
	float		fX								=	0;
	float		fY								=	0;
	int			nCnt							=	nPairCnt;
	
	//
	//LAPACK inputs
	//
	//Note that lapack always takes colum ordered matrix
	integer  rhs								=	1;
	integer  info;
	integer  pvt_1[MAX_TPS_ORDER]				=	{0};    //max order is 15
	integer  pvt_2[MAX_TPS_ORDER]				=	{0};	//	
	integer  n									=	nOrder;
	integer  lda								=	nOrder;
	const double PI								=   3.14159;

	doublereal a[MAX_TPS_ORDER*MAX_TPS_ORDER]	=	{0}; //used to solve for x
	doublereal b[MAX_TPS_ORDER]					=	{0};

	doublereal	A[MAX_TPS_ORDER*MAX_TPS_ORDER]	=	{0}; //used to solve for y	
	doublereal	B[MAX_TPS_ORDER]				=	{0};
	double		dDist							=	0;
	int			i								=	0;
	
	nCnt		=	nCnt>nMax?nMax:nCnt;
	nOrder		=	nCnt+3;
	n			=	nOrder;
	lda			=	nOrder;
	//prepare the coefficients

	//first three rows
	//a[0][0]-a[2][2] are zeros
	for(nCol=3;nCol<nOrder;nCol++) 
	{
		//first row
		a[nCol*nOrder  ]= 1;
		a[nCol*nOrder+1]= pMinT[pPair[nCol-3].m_nTIdx].m_nX;
		a[nCol*nOrder+2]= pMinT[pPair[nCol-3].m_nTIdx].m_nY;
	}

	//for rest of the rows
	//
	for(nRow = 3; nRow< nOrder;nRow++)
	{
		a[nRow]				=	1;					//first col
		a[nOrder+nRow]		=	pMinT[pPair[nRow-3].m_nTIdx].m_nX;//second col
		a[2*nOrder+nRow]	=	pMinT[pPair[nRow-3].m_nTIdx].m_nY;//second col
		//
		//LHS
		//
		for(idx=0,nCol=3;nCol<nOrder;nCol++,idx++)
		{
			
			dDist				=	((pMinT[pPair[idx].m_nTIdx].m_nX-pMinT[pPair[nRow-3].m_nTIdx].m_nX)*
									(pMinT[pPair[idx].m_nTIdx].m_nX-pMinT[pPair[nRow-3].m_nTIdx].m_nX))+

									((pMinT[pPair[idx].m_nTIdx].m_nY-pMinT[pPair[nRow-3].m_nTIdx].m_nY)*
									(pMinT[pPair[idx].m_nTIdx].m_nY-pMinT[pPair[nRow-3].m_nTIdx].m_nY));
			a[nCol*nOrder+nRow] =	dDist*log(dDist+1e-5);
		}
		//
		//RHS
		//
		b[nRow]				=	pMinR[pPair[nRow-3].m_nRIdx].m_nX;
		B[nRow]				=	pMinR[pPair[nRow-3].m_nRIdx].m_nY;
	}
	//copy a to A
	memcpy(A,a,MAX_TPS_ORDER*MAX_TPS_ORDER*sizeof(doublereal));

	//
	//Solve the equation
	//
	dgesv_(&nOrder, &rhs, a, &lda, pvt_1, b, &lda, &info); //get X Txn
	if(info) return UTL_ERR_INTERNAL_ERROR;

	dgesv_(&nOrder, &rhs, A, &lda, pvt_2, B, &lda, &info); //get Y Txn
	if(info) return UTL_ERR_INTERNAL_ERROR;

	//
	//Check for singular solutions in final version!
	//

	//record the transformation parameters
	pTx->m_nCnt = nOrder;

	for(i=0;i<nOrder;i++)
	{
		pTx->m_dTxX[i] = b[i];
		pTx->m_dTxY[i] = B[i];
	}

	//record the control points
	for(i=0;i<nCnt;i++)
	{
		pTx->m_minCtrl[i] = pMinT[pPair[i].m_nTIdx];
	}
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLTransformMinutiaSet
-------------------------------------------------------------------------------------------------------------*/
int UTLTransformMinutiaSet(Minutia pMin[],const int nCnt,const TxParams* pTx)
{
	/*sanity check*/
	Minutia tmp[MAX_MINUTIAE];
	int		i					=	0;
	int		j					=	0;
	double  dDist				=	0;
	double  dX					=	0;
	double	dY					=	0;

	if(pMin == NULL || nCnt < 1 || pTx == NULL)
		return UTL_ERR_INVALID_PARAMETER;
	
	for(i=0;i<nCnt;i++)
	{
		dX	=	pTx->m_dTxX[0]+pTx->m_dTxX[1]*pMin[i].m_nX + pTx->m_dTxX[2]*pMin[i].m_nY;
		dY  =   pTx->m_dTxY[0]+pTx->m_dTxY[1]*pMin[i].m_nX + pTx->m_dTxY[2]*pMin[i].m_nY;
		for(j=3;j< pTx->m_nCnt;j++)
		{
			dDist	=	(pMin[i].m_nX-pTx->m_minCtrl[j-3].m_nX)*(pMin[i].m_nX-pTx->m_minCtrl[j-3].m_nX)+
						(pMin[i].m_nY-pTx->m_minCtrl[j-3].m_nY)*(pMin[i].m_nY-pTx->m_minCtrl[j-3].m_nY);
			dDist	=	dDist*log(dDist+1e-5);
			
			dX		+=	pTx->m_dTxX[j]*dDist;
			dY		+=	pTx->m_dTxY[j]*dDist;
		}
		tmp[i].m_nX		=	(int)dX;
		tmp[i].m_nY		=	(int)dY;
		tmp[i].m_nTheta	=	0; //do something about this later
	}
	//copy it back
	for(i=0;i<nCnt;i++)
	{
		pMin[i] = tmp[i];
	}
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
UTLGetMatchScore
-------------------------------------------------------------------------------------------------------------*/
int UTLGetMatchScore(const Minutia pRMin[],const int nRCnt,const Minutia pTMin[], const int nTCnt,double* pdScore)
{
	int		i		=	0;
	int		j		=	0;
	double	dR		=	0;
	int		nRFlags[MAX_MINUTIAE]	=	{0};
	int		nTFlags[MAX_MINUTIAE]	=	{0};
	double	dScore	=	0;
	int		nULX	=	1e4;	/*upper left x,y*/
	int		nULY	=	1e4;
	int		nBRX	=	0;	/*bottom right x,y*/
	int		nBRY	=	0;
	int		nX		=	0;
	int		nY		=	0;
	int		nRMatches=	0;
	int		nTMatches=	0;
	int		nRTotal	=	0;
	int		nTTotal	=	0;

	*pdScore = 0;

	for(i=0;i<nRCnt;i++)
	{
		for(j=0;j<nTCnt;j++)
		{
			if(nTFlags[j] || nRFlags[i]) /*pair already matched*/
				continue;


			nX		=	abs(pRMin[i].m_nX-pTMin[j].m_nX);
			nY		=	abs(pRMin[i].m_nY-pTMin[j].m_nY);
			if(nX <= XTOL && nY <= YTOL)

			/*dR		=	(pTMin[j].m_nX-pRMin[i].m_nX)*(pTMin[j].m_nX-pRMin[i].m_nX)+
						(pTMin[j].m_nY-pRMin[i].m_nY)*(pTMin[j].m_nY-pRMin[i].m_nY);

			if(dR <= RTOL*RTOL)*/
			{
				nTFlags[j]	=	1;
				nRFlags[i]	=	1;
				dScore++;
				/*update bounding box*/
				if(nULX > pRMin[i].m_nX)
					nULX	=	pRMin[i].m_nX;
				if(nULY > pRMin[i].m_nY)
					nULY	=	pRMin[i].m_nY;	

				if(nBRX < pRMin[i].m_nX)
					nBRX	=	pRMin[i].m_nX;
				if(nBRY < pRMin[i].m_nY)
					nBRY	=	pRMin[i].m_nY;
			}
		}
	}
	
	/*check matches and mismatches*/
	for(i=0;i<nRCnt;i++)
	{
		nX	= pRMin[i].m_nX;
		nY	= pRMin[i].m_nY;
		if(nX <= nBRX && nX >= nULX && nY <= nBRY && nY >= nULY)
		{
			if(nRFlags[i])
				nRMatches+=4;
			else
				nRMatches-=2;
			nRTotal++;
		}
		
	}
	nRMatches = (nRMatches>0)?nRMatches/4:0;

	for(i=0;i<nTCnt;i++)
	{
		nX	= pTMin[i].m_nX;
		nY	= pTMin[i].m_nY;
		if(nX <= nBRX && nX >= nULX && nY <= nBRY && nY >= nULY)
		{
			if(nTFlags[i])
				nTMatches+=4;
			else
				nTMatches-=2;
			nTTotal++;
		}
		
	}
	nTMatches = (nTMatches>0)?nTMatches/4:0;
	printf("(%d,%d),(%d,%d)\n",nRMatches,nTMatches,nRTotal,nTTotal);
	//*pdScore = (double)(nTMatches*nRMatches)/(nRCnt*nTCnt)*100;//(nRCnt*nTCnt);
	*pdScore = (dScore*dScore)/(nRCnt*nTCnt)*100;
	return UTL_ERR_NO_ERROR;
}

/*-------------------------------------------------------------------------------------------------------------
MatchMinutiaSets
-------------------------------------------------------------------------------------------------------------*/
int UTLGetMatchScoreRadial(const RMinutia pRMin[],const int nRCnt,const RMinutia pTMin[],const int nTCnt,double* pdScore)
{
	int i			=	0;
	int j			=	0;
	double  dR		=	0;
	double	dTheta	=	0;
	double	dDelta	=	0;
	int		nAbs	=	0;
	int		nTFlags[MAX_MINUTIAE]	={0};
	int		nRFlags[MAX_MINUTIAE]	={0};
	double  dScore	=	0;

	/*sanity check*/
	if(pRMin == NULL || nRCnt < 1 || pTMin == NULL || nTCnt < 1 || pdScore == NULL)
		return UTL_ERR_INVALID_PARAMETER;

	for(i=0;i<nRCnt;i++)
	{
		if(nRFlags[i])
		{
			DBGPrintf("===>%d continuining..\n",i);
			continue;
		}

		for(j=0;j<nTCnt;j++)
		{
			if(nTFlags[j])continue; //already matched

			/*check for match*/
			dR		=	fabs(log10((pRMin[i].m_dR+1e-5)/(pTMin[j].m_dR+1e-5)));
			
			nAbs	=	abs(pRMin[i].m_nTheta - pTMin[j].m_nTheta);
			if((360-nAbs) < nAbs)
				dTheta = (360-nAbs);
			else
				dTheta = nAbs;

			nAbs	=	abs(pRMin[i].m_nDelta - pTMin[j].m_nDelta);


			if(dR <= DR && dTheta <= DTHETA && !nRFlags[i]) //match found
			{
				nRFlags[i]	=	1;
				nTFlags[j]	=	1;
				//add vote to pair matrix
				dScore++;
			}
		}
	}
	*pdScore = dScore;
	return UTL_ERR_NO_ERROR;

}

/*-------------------------------------------------------------------------------------------------------------
MatchMinutiaSets
-------------------------------------------------------------------------------------------------------------*/
int MatchMinutiaSets(const int pnAX[], const int pnAY[], const int pnATheta[], int nACnt,
					 const int pnBX[], const int pnBY[], const int pnBTheta[], int nBCnt,
					 double *pdScore)
{

	const int	PAIR_CNT		=	TPS_PAIRS;
	Minutia		minR[MAX_MINUTIAE];
	Minutia		minT[MAX_MINUTIAE];
	Minutia		minTmp[MAX_MINUTIAE];

	Pair		pairs[MAX_MINUTIAE*MAX_MINUTIAE];

	TxParams	tx;

	int		nRCnt;
	int		nTCnt;
	int		nPairCnt	=	0;
	int		nIter		=	0;
	int		i			=	0;
	int		j			=	0;
	int		nOldCnt		=	0;

	if(nACnt < nBCnt) //shorter one is the reference
	{
		UTLCreateMinutiaList(pnAX,pnAY,pnATheta,nACnt,minR);
		UTLCreateMinutiaList(pnBX,pnBY,pnBTheta,nBCnt,minT);
		nRCnt = nACnt;
		nTCnt = nBCnt;
	}
	else
	{
		UTLCreateMinutiaList(pnAX,pnAY,pnATheta,nACnt,minT);
		UTLCreateMinutiaList(pnBX,pnBY,pnBTheta,nBCnt,minR);
		nRCnt = nBCnt;
		nTCnt = nACnt;
	}

	//check for minimum count of minutia
	if(nRCnt < PAIR_CNT || nTCnt < PAIR_CNT)
	{
		*pdScore = 0;
		DBGPrintf("MatchMinutiaSets: Insufficient points(%d,%d)\n",nRCnt,nTCnt);
		return UTL_ERR_INSUFFICIENT_POINTS;
	}
	//reset pairs
	memset(pairs,0,MAX_MINUTIAE*MAX_MINUTIAE*sizeof(Pair));


	/*get the initial matches*/
	UTLGetInitialPairs(minR,nRCnt,minT,nTCnt,pairs);

	nOldCnt	=	-1;
	nPairCnt  = PAIR_CNT;
	printf("-----------------------------\n");
	printf("Counts are are %d,%d\n",nRCnt,nTCnt);
	/*perform iterative matching*/
	while(nOldCnt != nPairCnt && nIter < 3)
	{
		nOldCnt = nPairCnt;
		/*find the transformation parameters*/
		UTLGetTransformationParams(minR,minT,pairs,nPairCnt,&tx);
		for(i=0;i<nTCnt;i++)
			minTmp[i]=minT[i];
		UTLTransformMinutiaSet(minTmp,nTCnt,&tx);
		memset(pairs,0,MAX_MINUTIAE*MAX_MINUTIAE*sizeof(Pair));
		UTLFindMatches(minR,nRCnt,minTmp,nTCnt,pairs,&nPairCnt);
		UTLSortPairs(pairs,nRCnt*nTCnt);
		printf("Scores are %d,%d\n",nOldCnt,nPairCnt);
		nIter++;
	}
	/*get the matching score*/
	UTLGetMatchScore(minR,nRCnt,minTmp,nTCnt,pdScore);
	return UTL_ERR_NO_ERROR;
}
/*-------------------------------------------------------------------------------------------------------------
UTLFindMatches
-------------------------------------------------------------------------------------------------------------*/
int UTLFindMatches(Minutia pMinR[],int nRCnt,Minutia pMinT[],int nTCnt,Pair pPairs[],int* pnCnt)
{
	int i	=	0;
	int	j	=	0;
	int nRFlags[MAX_MINUTIAE]	=	{0};
	int nTFlags[MAX_MINUTIAE]	=	{0};
	int	nCnt					=	0;
	double dR					=	0;
	int	nDX						=	0;
	int nDY						=	0;

	*pnCnt =0;
	for(i=0;i<nRCnt;i++)
	{
		for(j=0;j<nTCnt;j++)
		{
			if(nTFlags[j] || nRFlags[i]) /*already matched*/
				continue;
			
			/*dR	=	(pMinR[i].m_nX-pMinT[j].m_nX)*(pMinR[i].m_nX-pMinT[j].m_nX) +
					(pMinR[i].m_nY-pMinT[j].m_nY)*(pMinR[i].m_nY-pMinT[j].m_nY);

			if(dR <= RTOL*RTOL)*/
			nDX		=	abs(pMinR[i].m_nX-pMinT[j].m_nX);
			nDY		=	abs(pMinR[i].m_nY-pMinT[j].m_nY);
			if(nDX <= XTOL && nDY <= YTOL)
			{
				nRFlags[i]	=	1;
				nTFlags[j]	=	1;
				nCnt++;
				pPairs[i*nTCnt+j].m_dScore +=	5;
				pPairs[i*nTCnt+j].m_nRIdx	=	i;
				pPairs[i*nTCnt+j].m_nTIdx	=	j;
			}
		}
	}
	*pnCnt = nCnt;
	return UTL_ERR_NO_ERROR;
}


