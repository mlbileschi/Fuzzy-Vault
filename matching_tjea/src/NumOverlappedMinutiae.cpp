#include "math.h"
#include "../inc/SFeature.h"

/*
NumOverlappedMinutiae()
Get the number of minutiae within the overlapped area of query and template prints.
Also return the width and height of combined fingerprints
*/
int NumOverlappedMinutiae(
                          const Minutiae   *arrM1, // For query fingerprint
                          const int        nMCnt1,
                          const int        refIdx1,
                          const Minutiae   *arrM2, // For template fignerprint
                          const int        nMCnt2,
                          const int        refIdx2,
                          const double     primary_od,
                          int              *pOLCnt1, // For query fingerprint
                          int              *pOLCnt2, // For template fingerprint
						  /*
                          int*             comb_w, // width of the combined fingerprint
                          int*             comb_h  // height of the combined fingerprint
						  */
						  MatchResultEx    *pResultEx // initialized by calling function
                          )
{
	int i;
    const Minutiae* arrM_A = NULL;
    Minutiae* arrM_B = NULL;
    int delta_x = 0, delta_y = 0;
    double cos_theta, sin_theta;
    struct {
        int minx;
        int maxx;
        int miny;
        int maxy;
    }comb_bb, a_bb, b_bb; // The bounding box of combined fingerprint
    
    // Initialize comb_bb
    a_bb.minx = b_bb.minx = comb_bb.minx = 10000;
    a_bb.miny = b_bb.miny = comb_bb.miny = 10000;
    a_bb.maxx = b_bb.maxx = comb_bb.maxx = -1;
    a_bb.maxy = b_bb.maxy = comb_bb.maxy = -1;
    
    *pOLCnt1 = *pOLCnt2 = -1;
	arrM_A = arrM1;
    arrM_B = (Minutiae*)malloc(sizeof(Minutiae) * nMCnt2);
    if(arrM_B == NULL) return -1;

    // Convert coordinates on sfv_arr2 into sfv_arr1 space
    cos_theta = cos(primary_od/180*PI);
    sin_theta = sin(primary_od/180*PI);
    delta_x = (arrM2[refIdx2].m_nX * cos_theta - arrM2[refIdx2].m_nY * sin_theta) - arrM1[refIdx1].m_nX;
    delta_y = (arrM2[refIdx2].m_nX * sin_theta + arrM2[refIdx2].m_nY * cos_theta) - arrM1[refIdx1].m_nY;
    
    for(i = 0; i < nMCnt2; i++){
		// find the bounding box for original (before convert coordinates) print
        if(arrM2[i].m_nX < b_bb.minx)
            b_bb.minx = arrM2[i].m_nX;
        else if(arrM2[i].m_nX > b_bb.maxx)
            b_bb.maxx = arrM2[i].m_nX;
        if(arrM2[i].m_nY < b_bb.miny)
            b_bb.miny = arrM2[i].m_nY;
        else if(arrM2[i].m_nY > b_bb.maxy)
            b_bb.maxy = arrM2[i].m_nY;
        arrM_B[i].m_nFeatureNo = arrM2[i].m_nFeatureNo;
        arrM_B[i].m_nX = (arrM2[i].m_nX * cos_theta - arrM2[i].m_nY * sin_theta) - delta_x;
        arrM_B[i].m_nY = (arrM2[i].m_nX * sin_theta + arrM2[i].m_nY * cos_theta) - delta_y;
        // find the bounding box
        if(arrM_B[i].m_nX < comb_bb.minx)
            comb_bb.minx = arrM_B[i].m_nX;
        else if(arrM_B[i].m_nX > comb_bb.maxx)
            comb_bb.maxx = arrM_B[i].m_nX;
        if(arrM_B[i].m_nY < comb_bb.miny)
            comb_bb.miny = arrM_B[i].m_nY;
        else if(arrM_B[i].m_nY > comb_bb.maxy)
            comb_bb.maxy = arrM_B[i].m_nY;
    }
    
    // find the bounding box
    for(i = 0; i < nMCnt1; i++){
        if(arrM_A[i].m_nX < comb_bb.minx)
            comb_bb.minx = arrM_A[i].m_nX;
        else if(arrM_A[i].m_nX > comb_bb.maxx)
            comb_bb.maxx = arrM_A[i].m_nX;
        if(arrM_A[i].m_nY < comb_bb.miny)
            comb_bb.miny = arrM_A[i].m_nY;
        else if(arrM_A[i].m_nY > comb_bb.maxy)
            comb_bb.maxy = arrM_A[i].m_nY;

		if(arrM1[i].m_nX < a_bb.minx)
            a_bb.minx = arrM1[i].m_nX;
        else if(arrM1[i].m_nX > a_bb.maxx)
            a_bb.maxx = arrM1[i].m_nX;
        if(arrM1[i].m_nY < a_bb.miny)
            a_bb.miny = arrM1[i].m_nY;
        else if(arrM1[i].m_nY > a_bb.maxy)
            a_bb.maxy = arrM1[i].m_nY;
    }
    /*
    *comb_w = comb_bb.maxx - comb_bb.minx + 1;
    *comb_h = comb_bb.maxy - comb_bb.miny + 1;
	*/
	pResultEx->Width_comb = comb_bb.maxx - comb_bb.minx + 1;
	pResultEx->Height_comb = comb_bb.maxy - comb_bb.miny + 1;
	pResultEx->Width_a = a_bb.maxx - a_bb.minx + 1;
	pResultEx->Height_a = a_bb.maxy - a_bb.miny + 1;
	pResultEx->Width_b = b_bb.maxx - b_bb.minx + 1;
	pResultEx->Height_b = b_bb.maxy - b_bb.miny + 1;


    *pOLCnt1 = NumMinutiaeInCH3(arrM_A, nMCnt1, arrM_B, nMCnt2, &b_bb.minx, &b_bb.miny, &b_bb.maxx, &b_bb.maxy);
    if(*pOLCnt1 < 0) 
    {
        free(arrM_B);
        return -1;
    }
    *pOLCnt2 = NumMinutiaeInCH3(arrM_B, nMCnt2, arrM_A, nMCnt1, &a_bb.minx, &a_bb.miny, &a_bb.maxx, &a_bb.maxy);
    if(*pOLCnt2 < 0)
    {
        free(arrM_B);
        return -1;
    }
    
    free(arrM_B);/*
	if((double)(a_bb.maxx-a_bb.minx)/(b_bb.maxx-b_bb.minx) > 1.5 ||
	   (double)(a_bb.maxx-a_bb.minx)/(b_bb.maxx-b_bb.minx) < 0.5 ||
	   (double)(a_bb.maxy-a_bb.miny)/(b_bb.maxy-b_bb.miny) > 1.5 ||
	   (double)(a_bb.maxy-a_bb.miny)/(b_bb.maxy-b_bb.miny) < 0.5)
	   return -2; //suspecious matching results.
	   */
    return 0;    //normal matching results.
}
#ifdef _noNEWSFV_
int NumOverlappedMinutiae(
                          const SFeature   *sfv_arr1, // For query fingerprint
                          const int        sfv_cnt1,
                          const int        refIdx1,
                          const SFeature   *sfv_arr2, // For template fignerprint
                          const int        sfv_cnt2,
                          const int        refIdx2,
                          const double     primary_od,
                          int              *pMinutiaeCnt1, // For query fingerprint
                          int              *pMinutiaeCnt2, // For template fingerprint
                          int*             comb_w, // width of the combined fingerprint
                          int*             comb_h  // height of the combined fingerprint
                          )
{
    int i;
    const SFeature* sfv_A = NULL;
    SFeature* sfv_B = NULL;
    int delta_x = 0, delta_y = 0;
    double cos_theta, sin_theta;
    struct {
        int minx;
        int maxx;
        int miny;
        int maxy;
    }comb_bb; // The bounding box of combined fingerprint
    
    // Initialize comb_bb
    comb_bb.minx = 10000;
    comb_bb.miny = 10000;
    comb_bb.maxx = -1;
    comb_bb.maxy = -1;
    
    *pMinutiaeCnt1 = *pMinutiaeCnt2 = -1;
    
    sfv_A = sfv_arr1;
    sfv_B = (SFeature*)malloc(sizeof(SFeature) * sfv_cnt2);
    if(sfv_B == NULL) return -1;
    
    // Convert coordinates on sfv_arr2 into sfv_arr1 space
    cos_theta = cos(primary_od/180*PI);
    sin_theta = sin(primary_od/180*PI);
    delta_x = (sfv_arr2[refIdx2].x * cos_theta - sfv_arr2[refIdx2].y * sin_theta) - sfv_arr1[refIdx1].x;
    delta_y = (sfv_arr2[refIdx2].x * sin_theta + sfv_arr2[refIdx2].y * cos_theta) - sfv_arr1[refIdx1].y;
    
    for(i = 0; i < sfv_cnt2; i++){
        sfv_B[i].FeatureNo = sfv_arr2[i].FeatureNo;
        sfv_B[i].x = (sfv_arr2[i].x * cos_theta - sfv_arr2[i].y * sin_theta) - delta_x;
        sfv_B[i].y = (sfv_arr2[i].x * sin_theta + sfv_arr2[i].y * cos_theta) - delta_y;
        // find the bounding box
        if(sfv_B[i].x < comb_bb.minx)
            comb_bb.minx = sfv_B[i].x;
        else if(sfv_B[i].x > comb_bb.maxx)
            comb_bb.maxx = sfv_B[i].x;
        if(sfv_B[i].y < comb_bb.miny)
            comb_bb.miny = sfv_B[i].y;
        else if(sfv_B[i].y > comb_bb.maxy)
            comb_bb.maxy = sfv_B[i].y;
    }
    
    // find the bounding box
    for(i = 0; i < sfv_cnt1; i++){
        if(sfv_A[i].x < comb_bb.minx)
            comb_bb.minx = sfv_A[i].x;
        else if(sfv_A[i].x > comb_bb.maxx)
            comb_bb.maxx = sfv_A[i].x;
        if(sfv_A[i].y < comb_bb.miny)
            comb_bb.miny = sfv_A[i].y;
        else if(sfv_A[i].y > comb_bb.maxy)
            comb_bb.maxy = sfv_A[i].y;
    }
    
    *comb_w = comb_bb.maxx - comb_bb.minx + 1;
    *comb_h = comb_bb.maxy - comb_bb.miny + 1;
    
    *pMinutiaeCnt1 = NumMinutiaeInCH2(sfv_A, sfv_cnt1, sfv_B, sfv_cnt2);
    if(*pMinutiaeCnt1 < 0) 
    {
        free(sfv_B);
        return -1;
    }
    *pMinutiaeCnt2 = NumMinutiaeInCH2(sfv_B, sfv_cnt2, sfv_A, sfv_cnt1);
    if(*pMinutiaeCnt2 < 0)
    {
        free(sfv_B);
        return -1;
    }
    
    free(sfv_B);
    return 0;
}
#endif /* noNEWSFV */
