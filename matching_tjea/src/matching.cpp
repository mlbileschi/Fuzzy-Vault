#pragma warning (disable : 4786)

#include <vector>
#include "../inc/common.h"
#include "../inc/matching.h"
#include "../inc/global_parameters.h"
#include "../inc/SFeature.h"
#include "../inc/logfile.h"

Params MatchingParams;

static void GetParams(){
	MatchingParams.dT_MaxSFVDist = Global_Parameters::M_T_MAXSFVDIST;
	MatchingParams.dT_MinSFVMatch = Global_Parameters::M_T_MINSFVMATCH;
	MatchingParams.nF_Hybrid = Global_Parameters::M_F_HYBRID;
	MatchingParams.nT_MinSFVReq = (int)Global_Parameters::M_T_MINSFVREQ;
	MatchingParams.nT_GoodWidth = (int)Global_Parameters::M_T_GOODWIDTH;
	MatchingParams.nT_GoodHeight = (int)Global_Parameters::M_T_GOODHEIGHT;
	MatchingParams.nT_TryNRef = (int)Global_Parameters::M_T_TRYNREF;
	MatchingParams.dT_F_MaxRadius = Global_Parameters::M_T_F_MAXRADIUS;
	MatchingParams.dT_P_MaxRadius = Global_Parameters::M_T_P_MAXRADIUS;
	MatchingParams.dT_F_MaxAngle = Global_Parameters::M_T_F_MAXANGLE;
	MatchingParams.dT_P_MaxAngle = Global_Parameters::M_T_P_MAXANGLE;
	MatchingParams.dT_F_MaxOrientation = Global_Parameters::M_T_F_MAXORIENTATION;
	MatchingParams.dT_P_MaxOrientation = Global_Parameters::M_T_P_MAXORIENTATION;
	MatchingParams.nF_BruteForce = (int)Global_Parameters::M_F_BRUTEFORCE;
}


int CUBS_MatchFeaturesInternal(ptrFPTemplate hRef, ptrFPTemplate hTst, 
							   MatchResult*	pMatchResult, MatchResultEx* pResultEx, std::vector<std::pair<int,int> > &matchingPairs)
{
//    MatchResult* pResult = NULL;
//	MatchResultEx* pResultEx = NULL;
    Minutiae* pMRef;
    Minutiae* pMTest;
    SFeature* pSFVRef = NULL;
    SFeature* pSFVTest = NULL;
	NeighborInfo NeighborListRef[CUBS_MAX_MINUTIAE];
	NeighborInfo NeighborListTest[CUBS_MAX_MINUTIAE];
	QUADSFVMAP RefQuadSFVList;       // SFV vectors for every quadrant of reference fingerprint
	QUADSFVMAP TestQuadSFVList;      // SFV vectors for every quadrant of test fingerprint


    int       nCntRef, nCntTest, i;
	int       nSFVCntRef;
	int       nSFVCntTest;
    int       ret = ERR_NO_ERROR;
    bool      MatchPerformed = false;
    double    similarity       = 0.0;
    int       pass = 1;

//!~	logfile_print("Beginning CUBS_MatchFeaturesInternal().\n");
    //Get parameters
    GetParams();

	if(hRef==NULL || hRef->nCnt<2 || hRef->nCnt>CUBS_MAX_MINUTIAE )
		return ERR_MISSING_MINUTIAE;
	pMRef   = (Minutiae*)hRef->marr;
    nCntRef = hRef->nCnt;
	for(i=0;i<nCntRef;i++)
		// matching algorithm asssumes there is an index here, 
		//not the type of minutia (as NIST minutia algorithm outputs)
		pMRef[i].m_nFeatureNo=i;  
	if(hTst==NULL || hTst->nCnt<2 || hTst->nCnt>CUBS_MAX_MINUTIAE )
		return ERR_MISSING_MINUTIAE;
	pMTest  = (Minutiae*)hTst->marr;
    nCntTest= hTst->nCnt;
	for(i=0;i<nCntTest;i++)
		// matching algorithm asssumes there is an index here, 
		//not the type of minutia (as NIST minutia algorithm outputs)
		pMTest[i].m_nFeatureNo=i;
   
    pass = 1;

    while(pass <= 1 /*MatchingParams.nT_Tries*/)
    {
		// pSFVRef has to be allocated inside Minutiae2SFeature
#ifdef FV_IDX
		// jea050305 start
		//ret = Minutiae2SFeature(nCntRef, pMRef, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntRef, &pSFVRef, NeighborListRef, RefQuadSFVList);
		if(nCntRef > 30)
			ret = Minutiae2SFeature(nCntRef, pMRef, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntRef, &pSFVRef, NeighborListRef, RefQuadSFVList);
		else if(nCntRef > 20 && nCntRef <= 30)
			ret = Minutiae2SFeature(nCntRef, pMRef, 7, 0, &nSFVCntRef, &pSFVRef, NeighborListRef, RefQuadSFVList);
		else
			ret = Minutiae2SFeature(nCntRef, pMRef, 10, 0, &nSFVCntRef, &pSFVRef, NeighborListRef, RefQuadSFVList);

		// jea050305 end
#else
		ret = Minutiae2SFeature(nCntRef, pMRef, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntRef, &pSFVRef, NeighborListRef);
#endif
//!~		logfile_print("CUBS_MatchFeaturesInternal(): extracted %d features from ref template\n", nSFVCntRef);
        if(ret != ERR_NO_ERROR){
            return ret;
        }

#ifdef FV_IDX
		// jea050305 start
        //ret = Minutiae2SFeature(nCntTest, pMTest, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntTest, &pSFVTest, NeighborListTest, TestQuadSFVList);
		if(nCntTest > 30)
			ret = Minutiae2SFeature(nCntTest, pMTest, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntTest, &pSFVTest, NeighborListTest, TestQuadSFVList);
		else if(nCntTest > 20 && nCntTest <=30)
			ret = Minutiae2SFeature(nCntTest, pMTest, 7, 0, &nSFVCntTest, &pSFVTest, NeighborListTest, TestQuadSFVList);
		else
			ret = Minutiae2SFeature(nCntTest, pMTest, 10, 0, &nSFVCntTest, &pSFVTest, NeighborListTest, TestQuadSFVList);
		// jea050305 end
#else
		ret = Minutiae2SFeature(nCntTest, pMTest, CUBS_NUM_OF_NEIGHBORS, 0, &nSFVCntTest, &pSFVTest, NeighborListTest);
#endif
//!~		logfile_print("CUBS_MatchFeaturesInternal(): extracted %d features from test template\n", nSFVCntTest);

        if(ret != ERR_NO_ERROR){
			free(pSFVRef);
            return ret;
        }
        
		memset(pResultEx, 0, sizeof(MatchResultEx));

        //----------------------------------------------------
        //Do some matching here
        //----------------------------------------------------
		pResultEx->NumOfMinutiae_a = nCntTest;
		pResultEx->NumOfMinutiae_b = nCntRef;
		ret = SFV_VerifyMatch(nCntTest, pMTest, NeighborListTest,
			                  nCntRef, pMRef, NeighborListRef,
			                  pSFVTest, nSFVCntTest,
#ifdef FV_IDX
							  TestQuadSFVList,
#endif
							  pSFVRef, nSFVCntRef,
#ifdef FV_IDX
							  RefQuadSFVList,
#endif
			                  &MatchPerformed, &similarity, pResultEx, matchingPairs);
        if(ret != ERR_NO_ERROR){
            free(pSFVRef);
            free(pSFVTest);
			free(pResultEx);
            return ret;
        }
        free(pSFVRef);
        free(pSFVTest);
        if(similarity > 0.0)
            break;
        pass++;
    }
    
    pMatchResult->MatchPerformed = MatchPerformed;
    pMatchResult->similarity = similarity;
		
//!~	logfile_print("CUBS_MatchFeaturesInternal(): similarity=%f\n", similarity);


    return ERR_NO_ERROR;
}
