/*!
@file   Minutiae2SFeature.cpp
@discussion
Provides implementations of the function Minutiae2SFeature
*/
#include <stdlib.h>
#include "../inc/matching.h"
#include "../inc/SFeature.h"
#include <math.h>
#include <algorithm>
#include <vector>
using namespace std;

// START------------ indexing sfv features -----------------
static bool IsValidQPair(QUADRANT q1, QUADRANT q2)
{
	if(q1 < QUADRANT_5)
		if((q2 - q1) >= 0)
			return ((q2 - q1) < 5);
		else return false;
	else
		if(((q2 - q1) + 8) >= 0)
			return (((q2 - q1) + 8) < 5);
		else return false;
}

static QUADRANT DecideQuadrant(double angle)
{
	// angle should be between 0 and 360
	if(angle < 45.0)      return QUADRANT_0;
	else if(angle < 90.0) return QUADRANT_1;
	else if(angle < 135.0)return QUADRANT_2;
	else if(angle < 180.0)return QUADRANT_3;
	else if(angle < 225.0)return QUADRANT_4;
	else if(angle < 275.0)return QUADRANT_5;
	else if(angle < 315.0)return QUADRANT_6;
	else                  return QUADRANT_7; // (angle < 360.0)
}

static QPAIR AddSFVtoQuadList(SFeature sfv, QUADSFVMAP& QuadSFVList)
{
	QUADRANT q1 = 0; // quadrant index for leg A
	QUADRANT q2 = 0; // quadrant index for leg B
	QUADRANT q1_alt = 0; // alternate quadrant index for leg A
	QUADRANT q2_alt = 0; // alternate quadrant index for leg B
	QPAIR    Q_Idx;

	q1 = DecideQuadrant(sfv.A.angle);
	q2 = DecideQuadrant(sfv.B.angle);
	Q_Idx = QPAIR(q1, q2);
	QuadSFVList[Q_Idx].push_back(sfv.FeatureNo); // index of secondary feature
	if((sfv.A.angle - q1) < 22.5)
		if(q1 == QUADRANT_0) q1_alt = QUADRANT_7;
		else q1_alt = q1 - 1;
	else
		if(q1 == QUADRANT_7) q1_alt = QUADRANT_0;
		else q1_alt = q1 + 1;
	if((sfv.B.angle - q2) < 22.5)
		if(q2 == QUADRANT_0) q2_alt = QUADRANT_7;
		else q2_alt = q2 - 1;
	else
		if(q2 == QUADRANT_7) q2_alt = QUADRANT_0;
		else q2_alt = q2 + 1;

	if(IsValidQPair(q1_alt, q2_alt))
	{
		QuadSFVList[QPAIR(q1_alt,q2_alt)].push_back(sfv.FeatureNo);
	}
	if(IsValidQPair(q1, q2_alt))
	{
		QuadSFVList[QPAIR(q1,q2_alt)].push_back(sfv.FeatureNo);
	}
	if(IsValidQPair(q1_alt, q2))
	{
		QuadSFVList[QPAIR(q1_alt,q2)].push_back(sfv.FeatureNo);
	}
	return Q_Idx;
}
// END-------------- indexing sfv features -----------------

static double dist_map[CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE];
// Indicates for every index, i, which column indexes, can be used to construct
// the SFV with i. (for mode 0)
// 
static unsigned char c_map[CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE];
// to indicate if a SFV already exist
// 0, if status is not set
// 1, if there is already a valid SFV
// 2, if a SFV is needed to be constructed (not sure if it is valid)
// 3, if it is not a valid SFV
static unsigned char status_map[CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE];

//-------------------------------------------------------------------
//MinutiaeCompare
//-------------------------------------------------------------------
static int MinutiaeCompare(const void* a, const void* b){
    if(((Minutiae*)a)->m_nY < ((Minutiae*)b)->m_nY) return -1;
    else if(((Minutiae*)a)->m_nY > ((Minutiae*)b)->m_nY) return 1;
    else{
        if(((Minutiae*)a)->m_nX < ((Minutiae*)b)->m_nX) return -1;
        else if(((Minutiae*)a)->m_nX > ((Minutiae*)b)->m_nX) return 1;
        else{
            if(((Minutiae*)a)->m_nTheta < ((Minutiae*)b)->m_nTheta){
                return -1;
            }else if(((Minutiae*)a)->m_nTheta > ((Minutiae*)b)->m_nTheta){
                return 1;
            }else return 0;
        }
    }
}

//-------------------------------------------------------------------
//cart_dist
//-------------------------------------------------------------------
static double cart_dist(CartPt A, CartPt B){
    return sqrt((double)((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y)));
}

//-------------------------------------------------------------------
//FindNearestNPoints
//-------------------------------------------------------------------
static void FindNearestNPoints(int       center_idx,
                               Minutiae* fv_arr,
                               int       pt_cnt,
                               int       N,
                               int       *out_idx_arr){
    int     i, j, k;
    CartPt  A, B;
    double  distance = 0.0;
    double *dist_cache = NULL;
    
    A.x = fv_arr[center_idx].m_nX;
    A.y = fv_arr[center_idx].m_nY;
    
    memset(out_idx_arr, -1, sizeof(int) * N);
    
    dist_cache = (double*)malloc(sizeof(double) * N);
    if(dist_cache == NULL){
        return;
    }
    for(i = 0; i < N; i++) dist_cache[i] = 65535.0;
    
    for(i = 0; i < pt_cnt; i++){
        if(i == center_idx) continue;
        B.x = fv_arr[i].m_nX;
        B.y = fv_arr[i].m_nY;
        distance = cart_dist(A,B);
        if(distance >= dist_cache[N - 1]){
            continue;
        }
        for(j = 0; j < N; j++){
            if(distance < dist_cache[j]){
                for(k = N - 1; k > j; k--){
                    out_idx_arr[k] = out_idx_arr[k - 1];
                    dist_cache[k] = dist_cache[k - 1];
                }
                out_idx_arr[j] = i;
                dist_cache[j] = distance;
                break;
            }
        }
    }
    free(dist_cache);
}

static unsigned char CHK_SFVStatus(int idx0, int idx1, int idx2)
{
	int root;
	int min_idx, max_idx;

	root = idx0;
	min_idx = idx1;
	if(idx2 < min_idx)
	{
		max_idx = min_idx;
		min_idx = idx2;
	}
	else
	{
		max_idx = idx2;
	}
	return status_map[root][min_idx][max_idx];
}

static unsigned char SET_SFVStatus(int idx0, int idx1, int idx2, unsigned char status)
{
	int root;
	int min_idx, max_idx;

	root = idx0;
	min_idx = idx1;
	if(idx2 < min_idx)
	{
		max_idx = min_idx;
		min_idx = idx2;
	}
	else
	{
		max_idx = idx2;
	}
	status_map[root][min_idx][max_idx] = status;
	return status_map[root][min_idx][max_idx];
}

static void GetNearestNIdx(int nMCnt, int i, int N)
{
	double min_dist = 65536.0;
	int    j = 0, n = N, idx = 0;
	
	for(n = N; n > 0; n--)
	{
		for(j = 0; j < nMCnt; j++)
		{
			if(i != j && c_map[i][j] == 0 && dist_map[i][j] < min_dist)
			{
				idx = j;
				min_dist = dist_map[i][j];
			}
		}
		c_map[i][idx] = 1;
		min_dist = 65536.0;
	}
}

static void CreateSFV(int i, int j, int k, Minutiae* arrM, int* pSFVCnt, SFeature* pSf)
{
	int x1, x2, y1, y2;
	int    v1x, v1y, v2x, v2y;
    double r = 0.0, theta = 0.0;

	pSf->RefNo = arrM[i].m_nFeatureNo;
	// this value should be assigned outside after removed the unwanted SFestures
	//pSf->FeatureNo = *pSFVCnt;
	(*pSFVCnt)++;
	pSf->x = arrM[i].m_nX;
	pSf->y = arrM[i].m_nY;
	pSf->orientation = arrM[i].m_nTheta;
	pSf->switched = 0;
	x1 = arrM[j].m_nX;
    y1 = arrM[j].m_nY;
    x2 = arrM[k].m_nX;
    y2 = arrM[k].m_nY;
    v1x = x1 - pSf->x;
    v1y = y1 - pSf->y;
	v2x = x2 - pSf->x;
	v2y = y2 - pSf->y;
	// Travers the angle COUNTER clockwise, the first met leg is A the other is B.
	if((v1x * v2y - v2x * v1y) >= 0){
		pSf->A.RefNo = arrM[j].m_nFeatureNo; // feature number of the minutia, should be the same as the index. if arrM is not sorted
		pSf->B.RefNo = arrM[k].m_nFeatureNo;
	}else{
		pSf->A.RefNo = arrM[k].m_nFeatureNo;
		pSf->B.RefNo = arrM[j].m_nFeatureNo;
	}
	
	cart2polar(arrM[pSf->A.RefNo].m_nX, arrM[pSf->A.RefNo].m_nY, pSf->x, pSf->y, &r, &theta);
	pSf->A.radius = r;
	// theta should be relative too (w.r.t. the orientation of center)
	//pSf->A.angle  = theta;
	pSf->A.angle  = theta - pSf->orientation;
	if(pSf->A.angle < 0) pSf->A.angle += 360.0;
	pSf->ra = arrM[pSf->A.RefNo].m_nTheta - arrM[i].m_nTheta;
	if(pSf->ra < 0)
		pSf->ra += 360;
	//pSf->A.orientation = arrM[pSf->A.RefNo].m_nTheta;
	pSf->A.orientation = pSf->ra; // should use relative info here
	pSf->A.RefNo = arrM[pSf->A.RefNo].m_nFeatureNo;

	
	cart2polar(arrM[pSf->B.RefNo].m_nX, arrM[pSf->B.RefNo].m_nY, pSf->x, pSf->y, &r, &theta);
	pSf->B.radius = r;
	// theta should be relative too (w.r.t. the orientation of center)
	//pSf->B.angle  = theta;
	pSf->B.angle  = theta - pSf->orientation;
	if(pSf->B.angle < 0) pSf->B.angle += 360.0;
	pSf->rb = arrM[pSf->B.RefNo].m_nTheta - arrM[i].m_nTheta;
	if(pSf->rb < 0)
		pSf->rb += 360;
	//pSf->B.orientation = arrM[pSf->B.RefNo].m_nTheta;
	pSf->B.orientation = pSf->rb;  // should use relative info here
	pSf->B.RefNo = arrM[pSf->B.RefNo].m_nFeatureNo;

	pSf->angle = pSf->B.angle - pSf->A.angle;
	if(pSf->angle < 0)
		pSf->angle += 360;
}

static bool pr(PolarPt p1, PolarPt p2)
{/*
	if(p1.angle == p2.angle)
	{
		if(p1.radius == p2.radius)
		{
			if(p1.orientation == p2.orientation)
			{
				return (p1.RefNo < p2.RefNo);
			}
			return (p1.orientation < p2.orientation);
		}
		return (p1.radius < p2.radius);
	}
	
	return (p1.angle < p2.angle);
	*/
	return (p1.radius < p2.radius);
}

static int GenerateAllPossibleSFV(int nMCnt, Minutiae* arrM, int mode, std::vector<SFeature> *pre_SFV, NeighborInfo* pNeighborList)
{
	int SFVCnt = 0;
	int i = 0, j = 0, k = 0;
	std::vector<SFeature>* pVec = pre_SFV;
	SFeature sf;

//	pVec = (std::vector<SFeature>*) pre_SFV;

//	pVec->clear();

	if(mode == 0) // create SFV from the central point and the closest n points
	{
		for(i = 0; i < nMCnt; i++)
		{
			PolarPt pt;
			pNeighborList[i].m = arrM[i];
			for(j = 0; j < nMCnt; j++)
			{
				for(k = 0; i != j && k < nMCnt; k++)
				{
					if(j != k && c_map[i][j] == 1 && c_map[i][k] == 1 && CHK_SFVStatus(i, j, k) == 0)
					{
						PtVector::iterator pt_it;
						SET_SFVStatus(i, j, k, 2);
						if(SFVCnt >= CUBS_MAX_FEATURES) return -1;// too many features
						CreateSFV(i, j, k, arrM, &SFVCnt, &sf);
						pt.RefNo = sf.A.RefNo;
						for(pt_it = pNeighborList[i].arrNeighbors.begin(); pt_it != pNeighborList[i].arrNeighbors.end(); pt_it++)
						{
							PolarPt pt2 = *pt_it;
							if(pt2.RefNo == pt.RefNo)
								break;
						}
						if(pt_it == pNeighborList[i].arrNeighbors.end()) // no identical item found
						{
							pt.orientation = sf.A.orientation;
							pt.radius      = sf.A.radius;
							pt.angle       = sf.A.angle;
							(pNeighborList[i]).arrNeighbors.push_back(pt);
						}
						pt.RefNo = sf.B.RefNo;
						for(pt_it = pNeighborList[i].arrNeighbors.begin(); pt_it != pNeighborList[i].arrNeighbors.end(); pt_it++)
						{
							PolarPt pt2 = *pt_it;
							if(pt2.RefNo == pt.RefNo)
								break;
						}
						if(pt_it == pNeighborList[i].arrNeighbors.end()) // no identical item found
						{
							pt.orientation = sf.B.orientation;
							pt.radius      = sf.B.radius;
							pt.angle       = sf.B.angle;
							(pNeighborList[i]).arrNeighbors.push_back(pt);
						}
						pVec->push_back(sf);
					}
				}//k
			}//j
			pNeighborList[i].nCnt = (pNeighborList[i].arrNeighbors).size();
			// Sort the arrNeighbors according to the angle.
			sort(pNeighborList[i].arrNeighbors.begin(),pNeighborList[i].arrNeighbors.end(), pr);
		}//i
		SFVCnt = pVec->size();
	}
	else if(mode == 1)
	{
	}

	return SFVCnt;
}

//-------------------------------------------------------------------
//Minutiae2SFeature
//Parameters:
//   [in] nMCnt    number of minutiae
//   [in] arrM     minutiae array
//   [in] N        number of closest points, if mode is 0
//                 radius to central point, if mode is 1
//   [in] mode     indicate how to generate SFV. (Usually, the SFV is 
//                 formed by the central point and its two closest points.
//                 That is mode=0 and N=2
//   [out] pnSVFCnt number of secondary features
//   [out] ppSFV    pointer of pointer of generated SFV array. Should allocated inside
//                 this function.
//-------------------------------------------------------------------
#ifdef FV_IDX
	int Minutiae2SFeature(int nMCnt, Minutiae* arrM, int N, int mode, int* pnSFVCnt, SFeature** ppSFV, NeighborInfo* pNeighborList, QUADSFVMAP& QuadSFVList){
#else
	int Minutiae2SFeature(int nMCnt, Minutiae* arrM, int N, int mode, int* pnSFVCnt, SFeature** ppSFV, NeighborInfo* pNeighborList){
#endif
	int		i = 0, j =0;
	CartPt  A, B;
	std::vector<SFeature> pre_SFV; // secondary features that prior to validating 
	std::vector<SFeature>::iterator it;
	
	if( nMCnt < 0 || arrM == NULL || nMCnt > CUBS_MAX_MINUTIAE || mode < 0 || mode > 1){
        return MATCH_ERR_INVALID_INPUT;
    }

	// Initialize output
	*pnSFVCnt = 0;
	*ppSFV = NULL;

	// Initialize dist_map (do not need)
	// ...
	// Initialize status_map
	memset(status_map, 0, sizeof(status_map));
	// Initialize c_map
	memset(c_map, 0, sizeof(c_map));
	
	// fill in dist_map
	for(i = 0; i < nMCnt; i++)
	{
		A.idx = arrM[i].m_nFeatureNo;
		A.x   = arrM[i].m_nX;
		A.y   = arrM[i].m_nY;
		for(j = i; j < nMCnt; j++)
		{
			B.idx = arrM[j].m_nFeatureNo;
			B.x   = arrM[j].m_nX;
			B.y   = arrM[j].m_nY;
			dist_map[i][j] = cart_dist(A,B);
			dist_map[j][i] = dist_map[i][j];
			if(mode == 1) // mark the SFV within radius N
			{
				if(i != j && dist_map[i][j] < (double)N)
				{
					c_map[i][j] = 1;
					c_map[j][i] = 1;
				}
			}
		} // j
		if(mode == 0)
		{
			GetNearestNIdx( nMCnt, i, N);
		}
	} // i

	// Create all possible SFV
	*pnSFVCnt = GenerateAllPossibleSFV(nMCnt, arrM, mode, &pre_SFV, pNeighborList);
	if(*pnSFVCnt <= 0) return ERR_TOO_MANY_FEATURES;

	// create SFV
	*ppSFV = (SFeature*)malloc(*pnSFVCnt * sizeof(SFeature));
	if(*ppSFV == NULL)
		return ERR_MEMORY_ERROR;
	*pnSFVCnt = 0;
	for(it = pre_SFV.begin(); it != pre_SFV.end(); it++)
	{
		// Only copy the valid SFV
		//SFeature *psfv;
		//psfv = it;
		if(it->angle < 10.0 || it->angle > 170.0) continue;
		// Now assign FeatureNo as array index
		it->FeatureNo = (*pnSFVCnt);
#ifdef FV_IDX
		psfv->Q_Idx = AddSFVtoQuadList(*psfv, QuadSFVList);
#endif
		//memcpy(*ppSFV+(*pnSFVCnt), it, sizeof(SFeature));
		*(*ppSFV+(*pnSFVCnt))=*it;
		(*pnSFVCnt)++;
	}
#ifdef FV_IDX
#ifdef CUBS_DEBUG
	// for debug information
	fprintf(stderr,"Total Secondary Features: %d\n", (*pnSFVCnt));
	for(QUADSFVMAP::iterator q_it = QuadSFVList.begin();q_it != QuadSFVList.end();q_it++)
	{
		fprintf(stderr,"Quadrant(%d,%d) = %d items\n", ((*q_it).first).first, ((*q_it).first).second, ((*q_it).second).size());
	}
#endif /* CUBS_DEBUG */
#ifdef QUADSTAT
	FILE* statFP = fopen("QuadStat.txt","a+");
	int max_cnt = 0;
	int min_cnt = 10000;
	int total_cnt = 0;
	int avg_cnt = 0;
	int binCnt  = 0;
	fprintf(statFP,"total %d\t", (*pnSFVCnt));
	for(QUADSFVMAP::iterator q_it = QuadSFVList.begin();q_it != QuadSFVList.end();q_it++)
	{
		if(((*q_it).second).size() > 0)
		{
			binCnt++;
			total_cnt += ((*q_it).second).size();
			if(((*q_it).second).size() > max_cnt) max_cnt = ((*q_it).second).size();
			if(((*q_it).second).size() < min_cnt) min_cnt = ((*q_it).second).size();
		}
	}
	fprintf(statFP,"avg\t%5.1f\tmin\t%d\tmax\t%d\n", (double)total_cnt/binCnt, min_cnt, max_cnt);
	fclose(statFP);
#endif /* QUADSTAT */
#endif /* FV_IDX */
	return ERR_NO_ERROR;
}
