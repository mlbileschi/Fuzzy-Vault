#include "../inc/SFeature.h"
#include "math.h"

// Copyright 2001, softSurfer (www.softsurfer.com)
// Modified by Tsai-Yang Jea on Dec. 4, 2003.
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.


// Assume that a class is already given for the object:
//    Point with coordinates {float x, y;}
//     (tjea:) We are using the CartPt {int x, int y} in SFeature.h
//===================================================================


// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm on Area of Triangles
inline int
isLeft( CartPt P0, CartPt P1, CartPt P2 )
{
    return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y - P0.y);
}

// isLeft2(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1 with 
//               considering the threshold thr
//            =0 for P2 on the line
//            <0 for P2 right of the line through P0 and P 1 with 
//               considering the threshold thr
// Author: Tsai-Yang Jea
inline int
isLeft2( CartPt P0, CartPt P1, CartPt P2, double thr )
{
	double det = (double)(P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y - P0.y);
	if(det >= 0.0)
		return (int) det;
	else
	{
		double r = sqrt((double)(P1.x - P0.x)*(P1.x - P0.x) + (P1.y - P0.y) * (P1.y - P0.y));
		double d = det/r;
		if(d + thr >= 0.0)
			return (int)(-1.0 * det);
		else
			return (int)det;
	}
}
//===================================================================


// chainHull_2D(): Andrew's monotone chain 2D convex hull algorithm
//     Input:  P[] = an array of 2D points 
//                   presorted by increasing x- and y-coordinates
//             n = the number of points in P[]
//     Output: H[] = an array of the convex hull vertices (max is n)
//     Return: the number of points in H[]
int
chainHull_2D( CartPt* P, int n, CartPt* H )
{
    // the output array H[] will be used as the stack
    int    bot=0, top=(-1);  // indices for bottom and top of the stack
    int    i;                // array scan index
    
    // Get the indices of points with min x-coord and min|max y-coord
    int minmin = 0, minmax;
    int xmin = P[0].x;
    for (i=1; i<n; i++)
        if (P[i].x != xmin) break;
        minmax = i-1;
        if (minmax == n-1) {       // degenerate case: all x-coords == xmin
            H[++top] = P[minmin];
            if (P[minmax].y != P[minmin].y) // a nontrivial segment
                H[++top] = P[minmax];
            H[++top] = P[minmin];           // add polygon endpoint
            return top+1;
        }
        
        // Get the indices of points with max x-coord and min|max y-coord
        int maxmin, maxmax = n-1;
        int xmax = P[n-1].x;
        for (i=n-2; i>=0; i--)
            if (P[i].x != xmax) break;
            maxmin = i+1;
            
            // Compute the lower hull on the stack H
            H[++top] = P[minmin];      // push minmin point onto stack
            i = minmax;
            while (++i <= maxmin)
            {
                // the lower line joins P[minmin] with P[maxmin]
                if (isLeft( P[minmin], P[maxmin], P[i]) >= 0 && i < maxmin)
                    continue;          // ignore P[i] above or on the lower line
                
                while (top > 0)        // there are at least 2 points on the stack
                {
                    // test if P[i] is left of the line at the stack top
                    if (isLeft( H[top-1], H[top], P[i]) > 0)
                        break;         // P[i] is a new hull vertex
                    else
                        top--;         // pop top point off stack
                }
                H[++top] = P[i];       // push P[i] onto stack
            }
            
            // Next, compute the upper hull on the stack H above the bottom hull
            if (maxmax != maxmin)      // if distinct xmax points
                H[++top] = P[maxmax];  // push maxmax point onto stack
            bot = top;                 // the bottom point of the upper hull stack
            i = maxmin;
            while (--i >= minmax)
            {
                // the upper line joins P[maxmax] with P[minmax]
                if (isLeft( P[maxmax], P[minmax], P[i]) >= 0 && i > minmax)
                    continue;          // ignore P[i] below or on the upper line
                
                while (top > bot)    // at least 2 points on the upper stack
                {
                    // test if P[i] is left of the line at the stack top
                    if (isLeft( H[top-1], H[top], P[i]) > 0)
                        break;         // P[i] is a new hull vertex
                    else
                        top--;         // pop top point off stack
                }
                H[++top] = P[i];       // push P[i] onto stack
            }
            if (minmax != minmin)
                H[++top] = P[minmin];  // push joining endpoint onto stack
            
            return top+1;
}

// To check if the point P in the convex hull H
// Return 1, if P is in or on H
// Return -1, if P is not in H
inline int IsInHull(const CartPt *H, const int n, CartPt P){
    int top = n - 1;
    while(top > 0){
        if(isLeft(H[top - 1], H[top], P) < 0){ // if P is right of the line H[top-1],H[top]
		            return -1;
        }
        top--;
    }
    if(isLeft(H[n - 1], H[0], P) < 0){
	    return -1;
    }
    return 1;
}

// To check if the point P in the convex hull H with considering threshold thr
// Return 1, if P is in or on H
// Return -1, if P is not in H
inline int IsInHull2(const CartPt *H, const int n, CartPt P, double thr){
    int top = n - 1;
    while(top > 0){
        if(isLeft2(H[top - 1], H[top], P, thr) < 0){ // if P is right of the line H[top-1],H[top]
		            return -1;
        }
        top--;
    }
    if(isLeft2(H[n - 1], H[0], P, thr) < 0){
	    return -1;
    }
    return 1;
}

static int PtCmp(const void* a, const void* b){
    CartPt* pa = (CartPt*) a;
    CartPt* pb = (CartPt*) b;
    if(pa->x == pb->x){
        return (pa->y - pb->y);
    }else{
        return (pa->x - pb->x);
    }
}

// Calaulate the number of minutiae in the Convex Hull of Matched minutiae set
// also returnds the bounding box of the convex hull
int NumMinutiaeInCH(const SFeature *sfv, const int SfvCnt, const int* MatchedList, const int MatchedCnt,
                    int* BB_minx, int* BB_miny, int* BB_maxx, int* BB_maxy){
    int i = 0;
    CartPt* PtList = NULL;
    CartPt* HullPts = NULL;
    CartPt* MatchedPts = NULL;
    int result = -1;
    int HullCnt = 0;
    
    if(MatchedCnt <= 0 || SfvCnt <= 0) return 0;
    
    *BB_minx = *BB_miny = 10000;
    *BB_maxx = *BB_maxy = -1;
    
    if(NULL == (PtList = (CartPt*)malloc(SfvCnt * sizeof(CartPt)))){
        return -1;
    }
    
    if(NULL == (HullPts = (CartPt*)malloc(SfvCnt * sizeof(CartPt)))){
        free(PtList);
        return -1;
    }
    
    if(NULL == (MatchedPts = (CartPt*)malloc(MatchedCnt * sizeof(CartPt)))){
        free(HullPts);
        free(PtList);
        return -1;
    }
    
    memset(PtList, -1, sizeof(CartPt) * SfvCnt);
    memset(HullPts, -1, sizeof(CartPt) * SfvCnt);
    memset(MatchedPts, -1, sizeof(CartPt) * MatchedCnt);
    
    for(i = 0; i < SfvCnt; i++){
        PtList[i].x = sfv[i].x;
        PtList[i].y = sfv[i].y;
        if(i < MatchedCnt){
            MatchedPts[i].x = sfv[MatchedList[i]].x;
            MatchedPts[i].y = sfv[MatchedList[i]].y;
        }
    }
    
    qsort(MatchedPts, MatchedCnt, sizeof(CartPt), PtCmp);
    
    HullCnt = chainHull_2D(MatchedPts, MatchedCnt, HullPts);
    
    // the last point on hull maybe the same as the first one.
    if(memcmp(HullPts + (HullCnt - 1), HullPts, sizeof(CartPt)) == 0 ){
        HullCnt--;
        memset(HullPts + HullCnt, -1, sizeof(CartPt));
    }
    
    result = 0;
    for(i = 0; i < SfvCnt; i++){
        if(IsInHull(HullPts, HullCnt, PtList[i]) > 0){
            result ++;
        }
    }
    
    for(i = 0; i < HullCnt; i++){
        if(HullPts[i].x < *BB_minx)
            *BB_minx = HullPts[i].x;
        else
            if(HullPts[i].x > *BB_maxx)
                *BB_maxx = HullPts[i].x;
            if(HullPts[i].y < *BB_miny)
                *BB_miny = HullPts[i].y;
            else
                if(HullPts[i].y > *BB_maxy)
                    *BB_maxy = HullPts[i].y;
    }
    
    free(PtList);
    free(HullPts);
    free(MatchedPts);
    return result;
}

// Calculate the number of minutiae of print A that are within the convex hull of print B
int NumMinutiaeInCH2(const SFeature *sfv_A, const int SfvCnt_A, 
                     const SFeature *sfv_B, const int SfvCnt_B
                     ){
    int i = 0;
    CartPt* PtListA = NULL;
    CartPt* HullPts = NULL; // Convex hull on B
    CartPt* PtListB = NULL;
    int result = -1;
    int HullCnt = 0;
    
    if(SfvCnt_B <= 0 || SfvCnt_A <= 0) return 0;
    
    if(NULL == (PtListA = (CartPt*)malloc(SfvCnt_A * sizeof(CartPt)))){
        return -1;
    }
    
    // Need to allocate one element more, since the first point would be the last one (duplicate)
    if(NULL == (HullPts = (CartPt*)malloc((SfvCnt_B + 1)* sizeof(CartPt)))){
        free(PtListA);
        return -1;
    }
    
    if(NULL == (PtListB = (CartPt*)malloc(SfvCnt_B * sizeof(CartPt)))){
        free(HullPts);
        free(PtListA);
        return -1;
    }
    
    memset(PtListA, -1, sizeof(CartPt) * SfvCnt_A);
    memset(HullPts, -1, sizeof(CartPt) * (SfvCnt_B + 1));
    memset(PtListB, -1, sizeof(CartPt) * SfvCnt_B);
    
    for(i = 0; i < SfvCnt_A; i++){
        PtListA[i].idx = sfv_A[i].FeatureNo;
        PtListA[i].x = sfv_A[i].x;
        PtListA[i].y = sfv_A[i].y;
    }
    
    for(i = 0; i < SfvCnt_B; i++){
        PtListB[i].idx = sfv_B[i].FeatureNo;
        PtListB[i].x = sfv_B[i].x;
        PtListB[i].y = sfv_B[i].y;
    }
    
    qsort(PtListB, SfvCnt_B, sizeof(CartPt), PtCmp);
    
    HullCnt = chainHull_2D(PtListB, SfvCnt_B, HullPts);
        
    // the last point on hull maybe the same as the first one.
    if(memcmp(HullPts + (HullCnt - 1), HullPts, sizeof(CartPt)) == 0 ){
        HullCnt--;
        memset(HullPts + HullCnt, -1, sizeof(CartPt));
    }

    /*
    fprintf(stderr,"Convex Hull points...\n");
    for(i = 0; i < HullCnt; i++){
        fprintf(stderr,"H(%d,%d,%d)\n", HullPts[i].idx, HullPts[i].x, HullPts[i].y);
    }
    */


    result = 0;
    //fprintf(stderr,"Points in overlapped area...\n");
    for(i = 0; i < SfvCnt_A; i++){
        if(IsInHull(HullPts, HullCnt, PtListA[i]) > 0){
            //fprintf(stderr,"%d\n", PtListA[i].idx);
            result ++;
        }
    }

#ifdef CUBS_DEBUG
    {
        CartPt* PtListO;  // points in overlapped area
        CartPt* OHullPts; // Overlapped hull points
        int     OHullCnt;
        int     PtListOCnt;
        int     j;
        PtListOCnt = result;
        if(PtListOCnt > 0){
            PtListO = (CartPt*) malloc(sizeof(CartPt) * (PtListOCnt + 1));
            // Need to allocate one element more, since the first point would be the last one (duplicate)
            OHullPts = (CartPt*) malloc(sizeof(CartPt) * (PtListOCnt + 1));
            for(i = 0, j = 0; i < SfvCnt_A; i++){
                if(IsInHull(HullPts, HullCnt, PtListA[i]) > 0){
                    PtListO[j] = PtListA[i];
                    j++;
                }
            }
            qsort(PtListO, PtListOCnt, sizeof(CartPt), PtCmp);
            OHullCnt = chainHull_2D(PtListO, PtListOCnt, OHullPts);
            fprintf(stderr,"Convex Hull of overlapped area.\n");
            for(i=0; i<OHullCnt; i++){
                fprintf(stderr,"%d\n",OHullPts[i].idx);
            }
            free(PtListO);
            free(OHullPts);
        }else{
            fprintf(stderr,"No point in overlapped area.\n");
        }
    }
#endif /* CUBS_DEBUG */
    
    free(PtListA);
    free(HullPts);
    free(PtListB);
    return result;
}

int NumMinutiaeInCH3(const Minutiae *arrM_A, const int Cnt_A, 
                     const Minutiae *arrM_B, const int Cnt_B,
					 int* BB_minx, int* BB_miny, int* BB_maxx, int* BB_maxy
                     ){
    int i = 0;
    CartPt* PtListA = NULL;
    CartPt* HullPts = NULL; // Convex hull on B
    CartPt* PtListB = NULL;
    int result = -1;
    int HullCnt = 0;
	*BB_minx = *BB_miny = 10000;
    *BB_maxx = *BB_maxy = -1;
    
    if(Cnt_B <= 0 || Cnt_A <= 0) return 0;
    
    if(NULL == (PtListA = (CartPt*)malloc(Cnt_A * sizeof(CartPt)))){
        return -1;
    }
    
    // Need to allocate one element more, since the first point would be the last one (duplicate)
    if(NULL == (HullPts = (CartPt*)malloc((Cnt_B + 1)* sizeof(CartPt)))){
        free(PtListA);
        return -1;
    }
    
    if(NULL == (PtListB = (CartPt*)malloc(Cnt_B * sizeof(CartPt)))){
        free(HullPts);
        free(PtListA);
        return -1;
    }
    
    memset(PtListA, -1, sizeof(CartPt) * Cnt_A);
    memset(HullPts, -1, sizeof(CartPt) * (Cnt_B + 1));
    memset(PtListB, -1, sizeof(CartPt) * Cnt_B);
    
    for(i = 0; i < Cnt_A; i++){
        PtListA[i].idx = arrM_A[i].m_nFeatureNo;
        PtListA[i].x = arrM_A[i].m_nX;
        PtListA[i].y = arrM_A[i].m_nY;
    }
    
    for(i = 0; i < Cnt_B; i++){
        PtListB[i].idx = arrM_B[i].m_nFeatureNo;
        PtListB[i].x = arrM_B[i].m_nX;
        PtListB[i].y = arrM_B[i].m_nY;
    }
    
    qsort(PtListB, Cnt_B, sizeof(CartPt), PtCmp);
    
    HullCnt = chainHull_2D(PtListB, Cnt_B, HullPts);
        
    // the last point on hull maybe the same as the first one.
    if(memcmp(HullPts + (HullCnt - 1), HullPts, sizeof(CartPt)) == 0 ){
        HullCnt--;
        memset(HullPts + HullCnt, -1, sizeof(CartPt));
    }

    result = 0;
    //fprintf(stderr,"Points in overlapped area...\n");
    for(i = 0; i < Cnt_A; i++){
        //if(IsInHull(HullPts, HullCnt, PtListA[i]) > 0){
		if(IsInHull2(HullPts, HullCnt, PtListA[i], 7.0) > 0){ // thr == 10 pixels
            //fprintf(stderr,"%d\n", PtListA[i].idx);
            result ++;
			if(PtListA[i].x >= *BB_maxx)
				*BB_maxx = PtListA[i].x;
			if(PtListA[i].x <= *BB_minx)
				*BB_minx = PtListA[i].x;
			if(PtListA[i].y >= *BB_maxy)
				*BB_maxy = PtListA[i].y;
			if(PtListA[i].y <= *BB_miny)
				*BB_miny = PtListA[i].y;
        }
    }
#ifdef CUBS_DEBUG
    {
        CartPt* PtListO;  // points in overlapped area
        CartPt* OHullPts; // Overlapped hull points
        int     OHullCnt;
        int     PtListOCnt;
        int     j;
        PtListOCnt = result;
        if(PtListOCnt > 0){
            PtListO = (CartPt*) malloc(sizeof(CartPt) * (PtListOCnt + 1));
            // Need to allocate one element more, since the first point would be the last one (duplicate)
            OHullPts = (CartPt*) malloc(sizeof(CartPt) * (PtListOCnt + 1));
            for(i = 0, j = 0; i < Cnt_A; i++){
                //if(IsInHull(HullPts, HullCnt, PtListA[i]) > 0){
				if(IsInHull2(HullPts, HullCnt, PtListA[i], 7.0) > 0){
                    PtListO[j] = PtListA[i];
                    j++;
                }
            }
            qsort(PtListO, PtListOCnt, sizeof(CartPt), PtCmp);
            OHullCnt = chainHull_2D(PtListO, PtListOCnt, OHullPts);
            fprintf(stderr,"Convex Hull of overlapped area.\n");
            for(i=0; i<OHullCnt; i++){
                fprintf(stderr,"%d\n",OHullPts[i].idx);
            }
            free(PtListO);
            free(OHullPts);
        }else{
            fprintf(stderr,"No point in overlapped area.\n");
        }
    }
#endif /* CUBS_DEBUG */
    
    free(PtListA);
    free(HullPts);
    free(PtListB);
    return result;
}
