#include <vector>
#include "../inc/matching.h"
#include "../inc/SFeature.h"
#include "../inc/logfile.h"

#include <iostream>
/*
* SFV_VerifyMatch
*/
int SFV_VerifyMatch(
					int nMCntTest, Minutiae* arrMTest, NeighborInfo* NeighborListTest,
					int nMCntRef, Minutiae* arrMRef, NeighborInfo* NeighborListRef,
					SFeature         *sfv_arr,// couldn't be const here, need do swith in sub-function
					const int        sfv_cnt,
					
#ifdef FV_IDX
					QUADSFVMAP       TestQuadSFVList,
#endif
					const SFeature   *Temp_sfv_arr, 
					const int        Temp_sfv_cnt,
#ifdef FV_IDX
					QUADSFVMAP       RefQuadSFVList,
#endif
					bool             *MatchingPerformed,
					double           *similarity,
					MatchResultEx    *pResultEx, // initialized by calling function
					std::vector<std::pair<int,int> > &matchingPairs
)
{
	int ret = ERR_NO_ERROR;
	//int MatchedListA[CUBS_MAX_FEATURES];
	//int MatchedListB[CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES];
	int *MatchedListA=(int*)malloc(sizeof(int)*CUBS_MAX_FEATURES);
	int *MatchedListB=(int*)malloc(sizeof(int)*CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES);
	int NumMatched = 0;
	//double ScoreList[CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES];
	double *ScoreList=(double*)malloc(sizeof(double)*CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES);
	int NumOfMatchedMinutiae = 0;
	double global_orientation = 0.0;
	int i;

	//int force_flag = 0;
	int force_flag = MatchingParams.nF_BruteForce;
	if(force_flag)
		fprintf(stderr,"Using Brute-Force matching.\n");
	
	*similarity = 0.0;
	*MatchingPerformed = false;
	
	logfile_print("Beginning SFV_VerifyMatch()\n");
	
	/* Initialize Matched lists */
//	memset(MatchedListA, -1, CUBS_MAX_FEATURES * sizeof(int));
//	memset(MatchedListB, -1, CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES * sizeof(int));
	for(i = 0; i < CUBS_MAX_FEATURES; i++){
		MatchedListA[i] = -1;
	}
	for(i = 0; i < CUBS_MAX_FEATURES * CUBS_MAX_CANDIDATES; i++){
		MatchedListB[i] = -1;
		ScoreList[i] = -1.0;
	}
	
#ifdef QUALITY_CHECK
	if(sfv_cnt < 10 || Temp_sfv_cnt < 10){
		free(MatchedListA);
		free(MatchedListB);
		free(ScoreList);
		return ERR_NO_ERROR;
	}
#endif /* QUALITY_CHECK */
	
	if(sfv_cnt <= MatchingParams.nT_MinSFVReq && Temp_sfv_cnt <= MatchingParams.nT_MinSFVReq) 
		if(MatchingParams.nF_Hybrid) force_flag = 1;
		
	if(!MatchingParams.nF_Hybrid || !force_flag){
#ifdef FV_IDX
		ret	= FeatureMatch(sfv_arr, sfv_cnt, TestQuadSFVList, Temp_sfv_arr, Temp_sfv_cnt, RefQuadSFVList, 
			CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES,
			&NumMatched, MatchedListA, MatchedListB, ScoreList);
#else
		ret	= FeatureMatch(sfv_arr, sfv_cnt, Temp_sfv_arr, Temp_sfv_cnt, CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES,
			&NumMatched, MatchedListA, MatchedListB, ScoreList);
#endif
		logfile_print("SFV_VerifyMatch(): got NumMatched=%d\n", NumMatched);
		if(ret != 0){
			free(MatchedListA);
			free(MatchedListB);
			free(ScoreList);
			return ret;
		}
#ifdef CUBS_DEBUG
		for(i=0; i < NumMatched; i++){
			int j;
			fprintf(stderr,"(%d %d %d)<-> ", sfv_arr[MatchedListA[i]].RefNo,sfv_arr[MatchedListA[i]].A.RefNo,sfv_arr[MatchedListA[i]].B.RefNo);
			for(j=0; j < CUBS_MAX_CANDIDATES && MatchedListB[i*CUBS_MAX_CANDIDATES + j] >= 0; j++){
				fprintf(stderr,"(%d %d %d) (%8.6f) ", Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].RefNo, Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].A.RefNo, Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].B.RefNo, ScoreList[i*CUBS_MAX_CANDIDATES + j]);
			}
			fprintf(stderr,"\n");
		}
#endif /* CUBS_DEBUG */
		
		global_orientation = PostProcess(sfv_arr, Temp_sfv_arr, CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES, &NumMatched, MatchedListA, MatchedListB, ScoreList);
		if(global_orientation <= 0.0) NumMatched = 0;
		logfile_print("SFV_VerifyMatch(): got global_orientation=%lf\n", global_orientation);
		
#ifdef CUBS_DEBUG
		fprintf(stderr, "OD = %lf\n", global_orientation);
		for(i=0; i < NumMatched; i++){
			int j;
			fprintf(stderr,"(%d %d %d)<-> ", sfv_arr[MatchedListA[i]].RefNo,sfv_arr[MatchedListA[i]].A.RefNo,sfv_arr[MatchedListA[i]].B.RefNo);
			for(j=0; j < CUBS_MAX_CANDIDATES && MatchedListB[i*CUBS_MAX_CANDIDATES + j] >= 0; j++){
				fprintf(stderr,"(%d %d %d) (%8.6f) ", Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].RefNo, Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].A.RefNo, Temp_sfv_arr[MatchedListB[i*CUBS_MAX_CANDIDATES + j]].B.RefNo, ScoreList[i*CUBS_MAX_CANDIDATES + j]);
			}
			fprintf(stderr,"\n");
		}
#endif /* CUBS_DEBUG */
		NumOfMatchedMinutiae = GetNumOfMatchedMinutiae(sfv_arr, sfv_cnt, NeighborListTest, nMCntTest, Temp_sfv_arr, Temp_sfv_cnt, NeighborListRef, nMCntRef, 
			CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES, NumMatched, MatchedListA, MatchedListB, ScoreList, global_orientation, pResultEx, MODE_REGULAR);
		if(NumOfMatchedMinutiae <= 0 && (sfv_cnt <= MatchingParams.nT_MinSFVReq || Temp_sfv_cnt <= MatchingParams.nT_MinSFVReq)){
			force_flag = 1;
		}
		logfile_print("SFV_VerifyMatch(): got NumOfMatchedMinutiae=%d\n", NumOfMatchedMinutiae);
		//}else{ // Forced
	} // regular matching

//std::cout << MatchedListA[i] << " " << MatchedListB[i*CUBS_MAX_CANDIDATES] << std::endl;

	if (force_flag) { // Forced
					  /*
					  NumOfMatchedMinutiae = GetNumOfMatchedMinutiae(sfv_arr, sfv_cnt, Temp_sfv_arr, Temp_sfv_cnt, 
					  CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES, NumMatched, MatchedListA, MatchedListB, ScoreList, global_orientation, pResultEx, MODE_FORCED);
		*/
		NumOfMatchedMinutiae = GetNumOfMatchedMinutiae(sfv_arr, sfv_cnt, NeighborListTest, nMCntTest, Temp_sfv_arr, Temp_sfv_cnt, NeighborListRef, nMCntRef, 
			CUBS_MAX_FEATURES, CUBS_MAX_CANDIDATES, NumMatched, MatchedListA, MatchedListB, ScoreList, global_orientation, pResultEx, MODE_FORCED);
	}
	
	if(NumOfMatchedMinutiae < 0) return MATCH_ERR_INTERNAL_ERROR;
		
#ifndef NEW_SCORING /* Use old scoring */
	{
		int OLMinutiaeA = 0, OLMinutiaeB = 0;
		//int comb_w = 0, comb_h = 0;
		
		if(ScoreList[0] < 0) ScoreList[0] = 0;
		
		if(NumOfMatchedMinutiae > 0){
		/*
		NumOverlappedMinutiae(arrMTest, nMCntTest, MatchedListA[0], 
		arrMRef, nMCntRef, MatchedListB[0],
		global_orientation, &OLMinutiaeA, &OLMinutiaeB, &comb_w, &comb_h);
			*/
			NumOverlappedMinutiae(arrMTest, nMCntTest, MatchedListA[0], 
				arrMRef, nMCntRef, MatchedListB[0],
				global_orientation, &OLMinutiaeA, &OLMinutiaeB, pResultEx);
			
			// added tjea 030305 --------------------
			if(OLMinutiaeA <= 0)
				OLMinutiaeA = NumOfMatchedMinutiae;
			if(OLMinutiaeB <= 0)
				OLMinutiaeB = NumOfMatchedMinutiae;
			// --------------------------------------
			pResultEx->NumOfMatchedMinutiae = NumOfMatchedMinutiae;
			pResultEx->NumOfMinutiae_a      = nMCntTest;
			pResultEx->NumOfMinutiae_b      = nMCntRef;
			pResultEx->NumOfOverlap_a       = OLMinutiaeA;
			pResultEx->NumOfOverlap_b       = OLMinutiaeB;
			pResultEx->DistScore            = ScoreList[0];
			pResultEx->PointScore           = (double)(NumOfMatchedMinutiae*NumOfMatchedMinutiae)/(OLMinutiaeA*OLMinutiaeB);
			
			// Add. 081804
			if(OLMinutiaeA < 5)
				OLMinutiaeA = 5;
			if(OLMinutiaeB < 5)
				OLMinutiaeB = 5;
			// end. 081804
			
			if(NumOfMatchedMinutiae < 7 && (pResultEx->Width_comb > MatchingParams.nT_GoodWidth || pResultEx->Height_comb > MatchingParams.nT_GoodHeight))
				*similarity = 0.0;
			else{
			/* try jea050305 start
			if(NumOfMatchedMinutiae >= 20 
			&& (double)NumOfMatchedMinutiae - (3.0/5.0*OLMinutiaeA) > 0.0 
			&& (double)NumOfMatchedMinutiae - (3.0/5.0*OLMinutiaeB) > 0.0)
			*similarity = 1.0 * ScoreList[0];
			else{
			*similarity = (double)(NumOfMatchedMinutiae*NumOfMatchedMinutiae)/(OLMinutiaeA*OLMinutiaeB)*ScoreList[0];
			if(*similarity > 1.0) *similarity = 1.0;
			}
				*/
				if(NumOfMatchedMinutiae >= 30 
					&& (double)NumOfMatchedMinutiae - (3.0/5.0*OLMinutiaeA) > 0.0 
					&& (double)NumOfMatchedMinutiae - (3.0/5.0*OLMinutiaeB) > 0.0)
					*similarity = 1.0 * ScoreList[0];
				else{
					*similarity = (double)(NumOfMatchedMinutiae*NumOfMatchedMinutiae)/(OLMinutiaeA*OLMinutiaeB)*ScoreList[0];
					if(*similarity > 1.0) *similarity = 1.0;
				}
				// jea050305 end
			}
		}
		else
			*similarity = 0.0;
	}
#else  /* defined NEW_SCORING */
	{
		int OLMinutiaeA = 0, OLMinutiaeB = 0;
		int comb_w = 0, comb_h = 0;
		
		if(ScoreList[0] < 0) ScoreList[0] = 0;
		
		if(NumOfMatchedMinutiae > 0){
			if(-2 == NumOverlappedMinutiae(arrMTest, nMCntTest, MatchedListA[0], 
				arrMRef, nMCntRef, MatchedListB[0],
				global_orientation, &OLMinutiaeA, &OLMinutiaeB, &comb_w, &comb_h))
			{
				// suspecious matching result
				// do something to re-verify them.
				// for now we just reduce the NumOfMatchedMinutiae into 1/2
				NumOfMatchedMinutiae /= 2;
				//tmp_num_match = NumOfMatchedMinutiae;
				//pResultEx->NoOfMatchedMinutiae = NumOfMatchedMinutiae;
			}
			
			pResultEx->NumOfMatchedMinutiae = NumOfMatchedMinutiae;
			pResultEx->NumOfMinutiae_a      = nMCntTest;
			pResultEx->NumOfMinutiae_b      = nMCntRef;
			pResultEx->NumOfOverlap_a       = OLMinutiaeA;
			pResultEx->NumOfOverlap_b       = OLMinutiaeB;
			pResultEx->DistScore            = ScoreList[0];
			pResultEx->PointScore           = (double)(NumOfMatchedMinutiae*NumOfMatchedMinutiae)/(OLMinutiaeA*OLMinutiaeB);
			
			if(OLMinutiaeA < NumOfMatchedMinutiae)
				OLMinutiaeA = NumOfMatchedMinutiae;
			if(OLMinutiaeB < NumOfMatchedMinutiae)
				OLMinutiaeB = NumOfMatchedMinutiae;
			
			if(NumOfMatchedMinutiae < 8 && (comb_w > MatchingParams.nT_GoodWidth || comb_h > MatchingParams.nT_GoodHeight))
				*similarity = 0.0;
			else{
				*similarity = (double)(NumOfMatchedMinutiae*NumOfMatchedMinutiae)/(OLMinutiaeA*OLMinutiaeB);
			}
		}
		else
			*similarity = 0.0;

#ifdef CUBS_DEBUG
		fprintf(stderr,"comb_w = %d, comb_h = %d, #match = %d, score = %8.6f \n", comb_w, comb_h, NumOfMatchedMinutiae, *similarity);
#endif /* CUBS_DEBUG */
	}
#endif /* NEW_SCORING */
	*MatchingPerformed = true;
	

// For fuzzy vault //
	for(i = 0; i < NumOfMatchedMinutiae; i++){
		double maxScore = -1000;
		int maxIndex = -1;
		for(int j=0; j<CUBS_MAX_CANDIDATES; j++){
			//std::cout << MatchedListA[i] << " " << MatchedListB[i*CUBS_MAX_CANDIDATES] << std::endl;
			//std::cout << " score: " << ScoreList[(i*CUBS_MAX_CANDIDATES) + j] << " ";
			if(ScoreList[(i*CUBS_MAX_CANDIDATES) + j] > maxScore){
				maxScore = ScoreList[(i*CUBS_MAX_CANDIDATES) + j];
				maxIndex = j;
			}
		}
		matchingPairs.push_back(std::make_pair(MatchedListA[i], MatchedListB[i*CUBS_MAX_CANDIDATES + maxIndex]));
	}
*similarity = NumOfMatchedMinutiae;
// End fuzzy vault //


	free(MatchedListA);
	free(MatchedListB);
	free(ScoreList);

	return ERR_NO_ERROR;
}
