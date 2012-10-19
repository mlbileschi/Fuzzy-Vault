#include <math.h>
#include <string.h>
//#include "MinCostAugmentation.h"
#include "../inc/SFeature.h"
#include "../inc/matching.h"
#include "../inc/logfile.h"
#include <limits>
#define oo        (std::numeric_limits<int>::max()) /* psudo infinitive number */

static double compute_sfv_dist(SFeature* pfv1, const SFeature* pfv2);

static double do_switch_match(const SFeature* pfv1, const SFeature* pfv2){
    double dist = -1.0;
    SFeature fv1;
    
    fv1.switched = 1;
    fv1.FeatureNo = pfv1->FeatureNo;
    fv1.orientation = pfv1->orientation;
    fv1.angle = pfv1->angle;
    fv1.x = pfv1->x;
    fv1.y = pfv1->y;
    // switch A and B
    memcpy(&(fv1.B), &(pfv1->A), sizeof(fv1.B));
    memcpy(&(fv1.A), &(pfv1->B), sizeof(fv1.A));
    fv1.ra = pfv1->rb;
    fv1.rb = pfv1->ra;
    
    dist = compute_sfv_dist(&fv1, pfv2);
    return dist;
}

static void do_switch(SFeature* pfv){
    PolarPt tmp_pt;
    int tmp_int;
    
    pfv->switched = 1;
    
    tmp_int = pfv->ra;
    pfv->ra = pfv->rb;
    pfv->rb = tmp_int;
    
    memcpy(&tmp_pt, &(pfv->A), sizeof(tmp_pt));
    memcpy(&(pfv->A), &(pfv->B), sizeof(pfv->A));
    memcpy(&(pfv->B), &tmp_pt, sizeof(pfv->B));
}

static double compute_sfv_dist(SFeature* pfv1, const SFeature* pfv2){
    double dist = 0.0;
    double ra1 = 0.0, rb1 = 0.0;
    double log_ra1 = 0.0, log_rb1 = 0.0;
    double ra2 = 0.0, rb2 = 0.0;
    double log_ra2 = 0.0, log_rb2 = 0.0;
    double d_ra = 0.0, d_rb = 0.0;
    double d_angle = 0.0;
    int    d_oa = 0, d_ob = 0;
#if 0    
    if((pfv1->switched != 1) && (
        (180.0 - pfv1->angle < 15.0 && 180.0 - pfv2->angle < 15.0) ||
        (pfv1->angle < 15.0 && pfv2->angle < 15.0))){
        dist = do_switch_match(pfv1, pfv2);
        if (dist > 0.0){
            // switch A and B and return
            do_switch(pfv1);
            if(dist > MatchingParams.dT_MaxSFVDist) return -1.0;
            return dist;
        }else{
            pfv1->switched = 0;
        }
    }
#endif
	{
		// variables for radius distance
		double tolerate_ratio = 0.1;
		double r_unit = 0.0;
		double r_dist = 0.0; // radius distance.
		// variables for angular distrance
		double max_r  = 210.0;
		double min_r  = 10.0;
		double lb     = 10.0; //5.0; jea043005
		double ub     = 20.0;
		double a_unit = 0.0;
		double a_dist = 0.0; // angular distance.
		// variables for orientation distance
		double o_unit = 0.0;
		double o_dist = 0.0; // orientation distance

		// compute the radius distance
		r_unit = pfv2->A.radius * tolerate_ratio;
		if(r_unit < 7.0) r_unit = 7.0; // jea043005
		r_dist = (ABSD(pfv1->A.radius - pfv2->A.radius))/r_unit;
		r_unit = pfv2->B.radius * tolerate_ratio;
		r_dist += (ABSD(pfv1->B.radius - pfv2->B.radius))/r_unit;
        r_dist /= 2.0;
		//if(r_dist > 2.0) return -1.0;
		if(r_dist > MatchingParams.dT_F_MaxRadius) return -1.0;

		// compute the angular and orientation distance
		//!!!! WATCHOUT !!!! We cannot compare the orientation of center minutiae
		if(pfv2->A.radius >= max_r)
			a_unit = lb;
		else if(pfv2->A.radius <= min_r)
            a_unit = ub;
		else
			a_unit = ub - ((pfv2->A.radius - min_r) * (ub - lb) / (max_r - min_r));
		a_dist = abs_angle_diff(pfv1->A.angle, pfv2->A.angle) / a_unit;
		o_unit = a_unit;
		o_dist = abs_angle_diff(pfv1->ra, pfv2->ra) / o_unit;
		if(pfv2->B.radius >= max_r)
			a_unit = lb;
		else if(pfv2->B.radius <= min_r)
            a_unit = ub;
		else
			a_unit = ub - ((pfv2->B.radius - min_r) * (ub - lb) / (max_r - min_r));
		a_dist += abs_angle_diff(pfv1->B.angle, pfv2->B.angle) / a_unit;
		o_unit = a_unit;
		o_dist += abs_angle_diff(pfv1->rb, pfv2->rb) / o_unit;
		a_dist /= 2.0;
		o_dist /= 2.0;
		//if(a_dist > 2.0) return -1.0;
		if(a_dist > MatchingParams.dT_F_MaxAngle) return -1.0;
		//if(o_dist > 2.0) return -1.0;
		if(o_dist > MatchingParams.dT_F_MaxOrientation) return -1.0;
		dist = r_dist + a_dist + o_dist;
	}

    if(dist > MatchingParams.dT_MaxSFVDist) return -1.0;
    return dist;
}

/* 
* If there is no match, this function will put the MatchList[0] = -1 and score[0] = -1.0
* Otherwise, there is at least a match.
*/
#ifdef FV_IDX
static void
FindMatches(SFeature           *fv,          // testing template
            const SFeature     *fv_list,     // reference template
            const unsigned int list_sz, 
			QUADSFVMAP         QuadSFVList,
            int                *MatchList,   // Indecies to fv_list
            double             *scores,      // The scores of the best match
            const int          max_candidates
            ){
    const SFeature* ptr = NULL;
	SFeature        SFea;
    double          dist = 0.0;
    int             j = 0;
    int             k = 0;
    int             numCand = 0;
    double          tmp_double;
    int             tmp_int;
    double          tmp_double2;
    int             tmp_int2;
    
    /* Initialize the result */
    MatchList[0] = -1;
    scores[0]    = -1.0;
	//---------------- SFV indexing -------------
	QUADSFVMAP::iterator q_it;  // iterator
	QPAIR                Q_Idx; // quadrant pair index
	SFVIDXVEC::iterator  i_it;  // index vector iterator
	SFVIDXVEC            idx_vec;
	Q_Idx = fv->Q_Idx;
	q_it = QuadSFVList.find(Q_Idx);
	if(q_it == QuadSFVList.end()) return;
	idx_vec = (*q_it).second;
	idx_vec.size();
	for (i_it = idx_vec.begin(); i_it != idx_vec.end(); i_it++){
		ptr = fv_list + (*i_it);
#ifdef CHEATING
        { 
            double angle = 0.0;
            if((180.0 - fv->angle < 15.0 && 180.0 - ptr->angle < 15.0) ||
                (fv->angle < 15.0 && ptr->angle < 15.0)){
                angle = fv->orientation - ptr->orientation;
                if(angle < 0.0) angle += 360.0;
                if(angle > CHEATING_ANGLE && angle < (360.0 - CHEATING_ANGLE)) continue;
            }else{
                angle = fv->orientation - ptr->orientation;
                if(angle < 0.0) angle += 360.0;
                if(angle > CHEATING_ANGLE && angle < (360.0 - CHEATING_ANGLE)) continue;
            }
        }
#endif /* CHEATING */
        dist = compute_sfv_dist(fv, ptr);
        if (dist < 0.0) continue;
        else{
            if(scores[0] < 0.0){
                scores[0]     = dist;
                MatchList[0]  = ptr->FeatureNo; // the feature index in SFV list
            }else{
                for(j = 0; j < max_candidates; j++){
                    if(dist < scores[j]){
                        tmp_double   = scores[j];
                        tmp_int      = MatchList[j];
                        scores[j]    = dist;
                        MatchList[j] = ptr->FeatureNo; // the feature index refer to the original list;
                        // Shift the array
                        for(k = j + 1; k < max_candidates; k++){
                            tmp_double2  = scores[k];
                            tmp_int2     = MatchList[k];
                            scores[k]    = tmp_double;
                            MatchList[k] = tmp_int;
                            tmp_double   = tmp_double2;
                            tmp_int      = tmp_int2;
                        }
                        break;
                    }
                }
            }
            numCand ++;
        }
    }
}
#else /* not use secondary feature indexing */
static void
FindMatches(SFeature           *fv, 
            const SFeature     *fv_list, 
            const unsigned int list_sz, 
            int                *MatchList,   // Indecies to fv_list
            double             *scores,      // The scores of the best match
            const int          max_candidates
            ){
    const SFeature* ptr = NULL;
	double          dist = 0.0;
    unsigned int    i = 0;
    int             j = 0;
    int             k = 0;
    int             numCand = 0;
    double          tmp_double;
    int             tmp_int;
    double          tmp_double2;
    int             tmp_int2;
    
    /* Initialize the result */
    MatchList[0] = -1;
    scores[0]    = -1.0;
	for ( i = 0, ptr = fv_list; i < list_sz; i++, ptr++){
#ifdef CHEATING
        { 
            double angle = 0.0;
            if((180.0 - fv->angle < 15.0 && 180.0 - ptr->angle < 15.0) ||
                (fv->angle < 15.0 && ptr->angle < 15.0)){
                angle = fv->orientation - ptr->orientation;
                if(angle < 0.0) angle += 360.0;
                if(angle > CHEATING_ANGLE && angle < (360.0 - CHEATING_ANGLE)) continue;
            }else{
                angle = fv->orientation - ptr->orientation;
                if(angle < 0.0) angle += 360.0;
                if(angle > CHEATING_ANGLE && angle < (360.0 - CHEATING_ANGLE)) continue;
            }
        }
#endif /* CHEATING */
        dist = compute_sfv_dist(fv, ptr);
        if (dist < 0.0) continue;
        else{
            if(scores[0] < 0.0){
                scores[0]     = dist;
                MatchList[0]  = i; // the feature index in SFV list
            }else{
                for(j = 0; j < max_candidates; j++){
                    if(dist < scores[j]){
                        tmp_double   = scores[j];
                        tmp_int      = MatchList[j];
                        scores[j]    = dist;
                        MatchList[j] = i; // the feature index refer to the original list;
                        // Shift the array
                        for(k = j + 1; k < max_candidates; k++){
                            tmp_double2  = scores[k];
                            tmp_int2     = MatchList[k];
                            scores[k]    = tmp_double;
                            MatchList[k] = tmp_int;
                            tmp_double   = tmp_double2;
                            tmp_int      = tmp_int2;
                        }
                        break;
                    }
                }
            }
            numCand ++;
        }
    }
}
#endif /* FV_IDX not use secondary feature indexing */
/*
* Returns ERR_NO_ERROR if success.
*/
int
FeatureMatch(SFeature       *list_a, 
             const int      list_a_sz,
#ifdef FV_IDX
			 QUADSFVMAP     TestQuadSFVList,
#endif
             const SFeature *list_b, 
             const int      list_b_sz,
#ifdef FV_IDX
			 QUADSFVMAP     RefQuadSFVList,
#endif
             const int      max_features,
             const int      max_candidates,
             int*           NumMatched,
             int*           MatchedListA,   // size = max_features, contains index of SFV arr.
             int*           MatchedListB,   // size = max_features * max_candidates
             double         *ScoreList      // size = max_features * max_candidates
             ){
    SFeature *ptr = NULL;
    unsigned int   i = 0, j = 0;
    double         *scores = NULL;
    
    /* Initialize result */
    *NumMatched = 0;
    if (MatchedListA == NULL ||
        MatchedListB == NULL){
        return MATCH_ERR_INVALID_INPUT; // error
    }
    memset(MatchedListA, -1, max_features * sizeof(int));
    memset(MatchedListB, -1, max_features * max_candidates * sizeof(int));
    for(i = 0; i < max_features * max_candidates; i++){
        ScoreList[i] = -1.0;
    }
    
    scores = (double*) malloc(sizeof(double) * max_candidates);
    if(NULL == scores){
        return ERR_MEMORY_ERROR; // error
    }
    
    for(i = 0, ptr = list_a; i < list_a_sz; i++, ptr++){
        for(j = 0 ; j < max_candidates; j++) scores[j] = -1.0;
#ifdef FV_IDX
        FindMatches(ptr, list_b, list_b_sz, RefQuadSFVList, &MatchedListB[*NumMatched * max_candidates], scores, max_candidates);
#else
		FindMatches(ptr, list_b, list_b_sz, &MatchedListB[*NumMatched * max_candidates], scores, max_candidates);
#endif
        if(MatchedListB[*NumMatched * max_candidates] < 0) continue;
        MatchedListA[*NumMatched] = i;
        if (ScoreList != NULL){
            for(j = 0; j < max_candidates && scores[j] >= 0.0; j++){
                ScoreList[*NumMatched * max_candidates + j] = scores[j];
            }
        }
        (*NumMatched) ++;
    }
    
    free (scores);

    return ERR_NO_ERROR;
}

/*
* Function : PostProcess
*/
double
PostProcess(const SFeature *list_a,
			const SFeature *list_b,
			const int      max_features,
            const int      max_candidates,
            int*           NumMatched,
            int*           MatchedListA,    // size = max_features
            int*           MatchedListB,    // size = max_features * max_candidates
            double         *ScoreList       // size = max_features * max_candidates
            ){
    double   global_orientation = 0;
    double   *angle_list = NULL;
    int i = 0, j = 0, k = 0;
    int best_idx, max_cnt = 0;
    int cnt = 0;
    int idx0, idx1;
    double angle = 0.0;
    int *NewMatchedListA;
    int *NewMatchedListB;
	double score_bins[NUM_OF_BINS]; // accumalate scores in every bin
    int bins[NUM_OF_BINS];  /* each element contains the hit counts of 10 degrees interval */
                            /* bins[0]  =>   0 ~ 9
                            bins[1]  =>  10 ~ 19
                            :               :
    bins[35] => 350 ~ 359 */
    unsigned char ex_bins[NUM_OF_BINS];  /* Exclusive bins to show which bin is increased by current sfv */
	int *marks;
    int angle_range = 0;
    //double* l_scores = NULL; // size = *NumMatched
	double* l_scores = NULL; // size = *NumMatched * max_candidate; make multiple choises
    
    if (*NumMatched <= 0) return 0;

	for(i = 0; i < NUM_OF_BINS; i++)
		score_bins[i] = 0.0;
    
    //l_scores = (double*)malloc(*NumMatched * sizeof(double));
	l_scores = (double*)malloc(*NumMatched * max_candidates * sizeof(double));
    if(l_scores == NULL){
        return -1.0;
    }
    
    marks = (int*) malloc(*NumMatched * max_candidates * sizeof(int));
    if (NULL == marks){
        free(l_scores);
        return -1.0;
    }
    angle_list = (double*) malloc(*NumMatched * sizeof(double) * max_candidates);
    if (NULL == angle_list){
        free(l_scores);
        free (marks);
        return -1.0;
    }
    
    memset(marks, -1, *NumMatched * max_candidates * sizeof(int));
    for(i = 0; i <  *NumMatched * max_candidates; i++){
        angle_list[i] = 400; // Invalid angle
    }
    
    memset(bins, 0, sizeof(int) * NUM_OF_BINS);
    
    angle_range = 360/NUM_OF_BINS;
    
    for (i = 0; i < *NumMatched; i++){
        memset(ex_bins, 0, sizeof(ex_bins));
        for(j = 0; j < max_candidates && MatchedListB[i * max_candidates + j] >= 0; j++){
            //angle = list_a[MatchedListA[i]].A.angle - list_b[MatchedListB[i * max_candidates + j]].A.angle;
            // the distortion of orientation is too huge
            angle = list_a[MatchedListA[i]].orientation - list_b[MatchedListB[i * max_candidates + j]].orientation;
            if(angle < 0) angle += 360.0;
            
            angle_list[i * max_candidates + j] = angle;
            marks[i * max_candidates + j] = (((int)angle)/angle_range%NUM_OF_BINS);
            if(ex_bins[marks[i * max_candidates + j]] == 0){
                bins[marks[i * max_candidates + j]] ++;
                ex_bins[marks[i * max_candidates + j]] ++;
				score_bins[marks[i * max_candidates + j]] += ScoreList[i * max_candidates + j];
            }
        }
    }
    
	// analyze the histogram (BIN), and find the peak with its neighbors
	double best_score = 10000.0;
    for (i = 0; i < NUM_OF_BINS; i++){
#ifdef CUBS_DEBUG
        fprintf(stderr,"BIN[%d]: %d\n", i, bins[i]);
#endif /* CUBS_DEBUG */
        if(i == 0){
            idx0 = NUM_OF_BINS - 1; // left neighbor
            idx1 = i + 1;           // right neighbor
        }else if(i == (NUM_OF_BINS - 1)){
            idx0 = i - 1;
            idx1 = 0;
        }else{
            idx0 = i - 1;
            idx1 = i + 1;
        }
        if(max_cnt < bins[i] + bins[idx0] + bins[idx1]){
            max_cnt = bins[i] + bins[idx0] + bins[idx1];
            best_idx = i;
			best_score = score_bins[i] + score_bins[idx0] + score_bins[idx1];
        }else if(max_cnt == bins[i] + bins[idx0] + bins[idx1]){
            if(best_score > score_bins[i] + score_bins[idx0] + score_bins[idx1]){
				//max_cnt = bins[i] + bins[idx0] + bins[idx1];
				best_idx = i;
				best_score = score_bins[i] + score_bins[idx0] + score_bins[idx1];
			}
		}
    }
    
	// Make sure there are at least two SFVs matched
	if(max_cnt < 2 && *NumMatched > 1 && score_bins[best_idx] > 2.0)
	{
		free(marks);
        free(angle_list);
        free(l_scores);
		return -2.0;
	}

    if (best_idx == 0){
        idx0 = NUM_OF_BINS - 1;
        idx1 = best_idx + 1;
    }else if (best_idx == (NUM_OF_BINS - 1)){
        idx0 = best_idx - 1;
        idx1 = 0;
    }else{
        idx0 = best_idx - 1;
        idx1 = best_idx + 1;
    }
    
    //max_cnt = bins[idx0] + bins[best_idx] + bins[idx1]; // why need?
    
    NewMatchedListA = (int*) malloc(max_cnt * sizeof(int));
    if (NULL == NewMatchedListA){
        free(marks);
        free(angle_list);
        free(l_scores);
        return -1.0;
    }
    memset(NewMatchedListA, -1, max_cnt * sizeof(int));
    NewMatchedListB = (int*) malloc(max_cnt * max_candidates * sizeof(int)); // make multiple choise
    if (NULL == NewMatchedListB){
        free(marks);
        free(angle_list);
        free(NewMatchedListA);
        free(l_scores);
        return -1.0;
    }
    memset(NewMatchedListB, -1, max_cnt * max_candidates * sizeof(int));

	// loop through all the matched results. 
    for(i = 0; i < *NumMatched * max_candidates; i++){
        if(marks[i] == best_idx || marks[i] == idx0 || marks[i] == idx1){
            if(angle_list[i] > 180.0)//make it between -180 and 180?
                global_orientation += angle_list[i] - 360.0;
            else
                global_orientation += angle_list[i];
        }
    }
    global_orientation = global_orientation / ((double)max_cnt);
    if(global_orientation < 0)// convert it back to [0 360) ??
        global_orientation += 360;
    
	double* l_angles = (double*)malloc(max_candidates * sizeof(double));
    cnt = 0;
    for(i = 0, j = 0; i < *NumMatched; i++){
        int curr_idx = -1;
		for(k = 0; k < max_candidates && MatchedListB[i * max_candidates + k] >= 0; k++){
			curr_idx = i * max_candidates + k;
            if(marks[curr_idx] == (int)best_idx || 
               marks[curr_idx] == (int)idx0 || 
               marks[curr_idx] == (int)idx1){
				if(NewMatchedListA[j]==-1){
					cnt++;
					NewMatchedListA[j] = MatchedListA[i]; // need only once
					memset(l_angles, 0, sizeof(double)*max_candidates);
					NewMatchedListB[j*max_candidates+0] = MatchedListB[curr_idx];
					l_scores[j*max_candidates+0] = ScoreList[curr_idx];
					l_angles[0] = angle_list[curr_idx];
					continue;
				}
				for(int m = 0; m < max_candidates; m++)
				{
					if(NewMatchedListB[j*max_candidates+m] < 0 ||
					  (angle_list[curr_idx] != 400.0 && 
					   ScoreList[curr_idx] < l_scores[j*max_candidates+m])
					){
						for(int n = max_candidates - 1; n > m; n--){
							if(NewMatchedListB[j*max_candidates+n - 1] >= 0){
								NewMatchedListB[j*max_candidates+n] = NewMatchedListB[j*max_candidates+(n-1)];
								l_scores[j*max_candidates+n] = l_scores[j*max_candidates+(n-1)];
								l_angles[n] = l_angles[(n-1)];
							}
						}//n
						NewMatchedListB[j*max_candidates+m] = MatchedListB[curr_idx];
						l_scores[j*max_candidates+m] = ScoreList[curr_idx];
						l_angles[m] = angle_list[curr_idx];
						break;
					}
				} // m
            }
        }//k
        if(j < max_cnt && NewMatchedListB[j*max_candidates+0] >= 0 && NewMatchedListA[j] >= 0)
            j++;
    }
	free(l_angles);
    
	// reset outputs
    memset(MatchedListA, -1, max_features * sizeof(int));
    memset(MatchedListB, -1, max_features * max_candidates * sizeof(int));
	//ScoreList is reset later.
	
    *NumMatched = 0;
		
	// create output
	if (cnt >= MatchingParams.dT_MinSFVMatch){
		for(i = 0; i < max_features * max_candidates; i++) ScoreList[i] = -1.0;
		memcpy(MatchedListA, NewMatchedListA, cnt * sizeof(int));
		memcpy(MatchedListB, NewMatchedListB, cnt * max_candidates * sizeof(int));
		memcpy(ScoreList, l_scores, cnt * max_candidates * sizeof(double));
		*NumMatched = cnt;
	}else{
		return -1.0; // no match here
	}
	
    free(marks);
    free(angle_list);
    free(NewMatchedListA);
    free(NewMatchedListB);
    free(l_scores);

    if(*NumMatched <= 0) return -1.0;
    else
        return global_orientation;
}


/*
*  Function : FindRefPt
*/
static int 
FindRefPt(
          const SFeature     *list_a,
          const int          cnt_a,
          const SFeature     *list_b,
          const int          cnt_b,
          const int          max_features,
          const int          max_candidates,
          const int          NumMatched,
          int*               MatchedListA,    // size = max_features
          int*               MatchedListB,    // size = max_features * max_candidates
          double             *ScoreList,      // size = max_features * max_candidates
          int*               best_idx_A,
          int*               best_idx_B
          ){
    int ret = ERR_NO_ERROR; // OK
    int i = 0, j = 0;
    double best_score = 10000.0;
    
    *best_idx_A = -1;
    *best_idx_B = -1;
    
    if(NumMatched <= 0) return ret;
    
    for(i = 0; i < NumMatched; i++){
        if(ScoreList[i * max_candidates] < 0) break;
        for(j = 0; j < max_candidates; j++){
            if(ScoreList[i * max_candidates +j] < 0) break;
            if(ScoreList[i * max_candidates +j] < best_score){
                *best_idx_A = MatchedListA[i];
                *best_idx_B = MatchedListB[i * max_candidates +j];
                best_score = ScoreList[i * max_candidates +j];
            }
        }
    }
    
    return ret;
}

/********************************************************************
Function Name: GetNumOfMatchedMinutiae
Return:  Number of matched minutiae. ( >= 0)
-1 (if error odccurs)
Note: We only count the minutiae in matched secondary feature
now.  We need modify this to get more accurate number in
the future.         
*********************************************************************/
struct {int hits; double accScore;}corr_map[CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE]; // correspondence map between list_a and list_b
struct {int idx_in_b; int supports;} best_for_a[CUBS_MAX_MINUTIAE];
struct {int hits; double accScore;}corr_map_leg[CUBS_MAX_MINUTIAE][CUBS_MAX_MINUTIAE]; // correspondence map between list_a and list_b
struct {int idx_in_b; int supports;} best_for_a_leg[CUBS_MAX_MINUTIAE];

int GetNumOfMatchedMinutiae(
							const SFeature   *list_a,
							const int        cnt_a,          // number of secondary features
							NeighborInfo*    NeighborListA,
                            const int        cntMinutiae_a,  // number of minutiae
                            const SFeature   *list_b,
							const int        cnt_b,          // number of secondary features
							NeighborInfo*    NeighborListB,
                            const int        cntMinutiae_b,  // number of minutiae
                            const int        max_features,
                            const int        max_candidates,
                            int              NumMatched,
                            int*             MatchedListA,    // size = max_features
                            int*             MatchedListB,    // size = max_features * max_candidates
                            double           *ScoreList,      // size = max_features * max_candidates
                            double&          primary_od,      // Primary orientation difference, if it is in forced mode, we need to fill the value in this function
							                                  // Otherwise, it is an input value
							MatchResultEx    *pResultEx,
                            int              mode             // Regular mode or Forced mode
                            ){
	// In this function the contents in MatchedListA and MatchedListB will be converted from
	// the indeices of secondary feature lists into the indecies of miutiae lists.
	logfile_print("GetNumOfMatchedMinutiae(): step 0\n");

    int NumOfMatchedMinutiae = -1; // Initial to error case.
	int i = 0, j = 0, flag = 0;
	int root_a = 0, root_b = 0, legA_a = 0, legA_b = 0, legB_a = 0, legB_b = 0;
	double score = -1.0;
	int    total_hits = 0;
	double best_score = (double)oo;
	int    best_idx = -1;

	NumOfMatchedMinutiae = -1;
	// Initialize corr_map
	// jea043005 --start--
	/*
	memset(corr_map, 0, sizeof(corr_map));
	memset(best_for_a, 0, sizeof(best_for_a));
	memset(corr_map_leg, 0, sizeof(corr_map_leg));
	memset(best_for_a_leg, 0, sizeof(best_for_a_leg));
	*/
	memset(corr_map, 0, sizeof(corr_map));
	memset(best_for_a, -1, sizeof(best_for_a));
	memset(corr_map_leg, 0, sizeof(corr_map_leg));
	memset(best_for_a_leg, -1, sizeof(best_for_a_leg));
	// jea043005 --end--
	
	logfile_print("GetNumOfMatchedMinutiae(): step 1\n");

	for(i = 0; i < NumMatched; i++)
	{
		for(j = 0; j < max_candidates && MatchedListB[i*max_candidates+j] > 0; j++)
		{
			root_a = list_a[MatchedListA[i]].RefNo;
			root_b = list_b[MatchedListB[i*max_candidates+j]].RefNo;
			corr_map[root_a][root_b].hits++;
			corr_map[root_a][root_b].accScore += ScoreList[i*max_candidates+j];
			if(corr_map[root_a][root_b].hits > best_for_a[root_a].supports)
			{
				best_for_a[root_a].supports = corr_map[root_a][root_b].hits;
				best_for_a[root_a].idx_in_b = root_b;
			}
			legA_a = list_a[MatchedListA[i]].A.RefNo;
			legA_b = list_b[MatchedListB[i*max_candidates+j]].A.RefNo;
			corr_map_leg[legA_a][legA_b].hits++;
			corr_map_leg[legA_a][legA_b].accScore += ScoreList[i*max_candidates+j];
			if(corr_map_leg[legA_a][legA_b].hits > best_for_a_leg[legA_a].supports)
			{
				best_for_a_leg[legA_a].supports = corr_map_leg[legA_a][legA_b].hits;
				best_for_a_leg[legA_a].idx_in_b = legA_b;
			}
			legB_a = list_a[MatchedListA[i]].B.RefNo;
			legB_b = list_b[MatchedListB[i*max_candidates+j]].B.RefNo;
			corr_map_leg[legB_a][legB_b].hits++;
			corr_map_leg[legB_a][legB_b].accScore += ScoreList[i*max_candidates+j];
			if(corr_map_leg[legB_a][legB_b].hits > best_for_a_leg[legB_a].supports)
			{
				best_for_a_leg[legB_a].supports = corr_map_leg[legB_a][legB_b].hits;
				best_for_a_leg[legB_a].idx_in_b = legB_b;
			}
		}//j
	}// i
	logfile_print("GetNumOfMatchedMinutiae(): step 2\n");

	// reset output
	memset(MatchedListA, -1, sizeof(int) * max_features);
    memset(MatchedListB, -1, sizeof(int) * max_features * max_candidates);
	memset(ScoreList, -1, sizeof(double) * max_features * max_candidates);
	ScoreList[0] = 0.0;
	score = 0.0;
	
//#define jea050405
#ifndef jea050405
	// remove duplicates
	for(i = 0; i < CUBS_MAX_MINUTIAE; i++)
	{
		for(j = i + 1; best_for_a[i].supports > 0 && j < CUBS_MAX_MINUTIAE; j++)
		{
			if(best_for_a[j].idx_in_b == best_for_a[i].idx_in_b)
			{
				if(best_for_a[i].supports < best_for_a[j].supports)
				{
					best_for_a[i].supports = 0;
					best_for_a[i].idx_in_b = -1;//jea043005
					break;
				}
				else if(best_for_a[i].supports > best_for_a[j].supports)
				{
					best_for_a[j].supports = 0;
					best_for_a[j].idx_in_b = -1;//jea043005
				}/* keep both if they have same supports */
			}
		}
	}
	for(i = 0; i < CUBS_MAX_MINUTIAE; i++)
	{
		for(j = i + 1; best_for_a_leg[i].supports > 0 && j < CUBS_MAX_MINUTIAE; j++)
		{
			if(best_for_a_leg[j].idx_in_b == best_for_a_leg[i].idx_in_b)
			{
				if(best_for_a_leg[i].supports < best_for_a_leg[j].supports)
				{
					best_for_a_leg[i].supports = 0;
					best_for_a_leg[i].idx_in_b = -1;//jea043005
					break;
				}
				else if(best_for_a_leg[i].supports > best_for_a_leg[j].supports)
				{
					best_for_a_leg[j].supports = 0;
					best_for_a_leg[j].idx_in_b = -1;//jea043005
				}/* keep both if they have same supports */
			}
		}
	}
	for(i = 0; i < CUBS_MAX_MINUTIAE; i++)
	{
		// jea043005 --start--
		if(best_for_a[i].supports > 0 && best_for_a_leg[i].supports > 0 && best_for_a[i].idx_in_b == best_for_a_leg[i].idx_in_b) 
		{
			int total_supports = best_for_a[i].supports + best_for_a_leg[i].supports;
			flag = 0;
			for(j = 0; j < CUBS_MAX_MINUTIAE; j++)
			{
				if(i != j && best_for_a[i].idx_in_b == best_for_a[j].idx_in_b
					&& best_for_a[i].idx_in_b == best_for_a_leg[j].idx_in_b
					&& best_for_a[j].supports + best_for_a_leg[j].supports > total_supports)
				{
					flag = 1;
					break;
				}
			}
			if(flag) continue;
			NumOfMatchedMinutiae++;
			MatchedListA[NumOfMatchedMinutiae] = i;
			MatchedListB[NumOfMatchedMinutiae * max_candidates] = best_for_a[i].idx_in_b;
			ScoreList[NumOfMatchedMinutiae * max_candidates] = (corr_map[i][best_for_a[i].idx_in_b].accScore+corr_map_leg[i][best_for_a[i].idx_in_b].accScore)/(corr_map[i][best_for_a[i].idx_in_b].hits+corr_map_leg[i][best_for_a[i].idx_in_b].hits);
			
			score += corr_map[i][best_for_a[i].idx_in_b].accScore;
			total_hits += corr_map[i][best_for_a[i].idx_in_b].hits;
			score += (corr_map_leg[i][best_for_a_leg[i].idx_in_b].accScore);
			total_hits += (corr_map_leg[i][best_for_a_leg[i].idx_in_b].hits);
		}
        /*
		if(best_for_a[i].idx_in_b == best_for_a_leg[i].idx_in_b 
			&& best_for_a[i].supports > 0
			&& best_for_a_leg[i].supports > 0)
		{
			flag = 0;
			for(j = 0; j < CUBS_MAX_MINUTIAE; j++)
			{
				if(i != j && best_for_a[i].idx_in_b == best_for_a[j].idx_in_b
					&& best_for_a[i].idx_in_b == best_for_a_leg[j].idx_in_b
					&& best_for_a[j].supports > 0)
					flag = 1;
			}
			for(j = 0; flag != 1 && j < CUBS_MAX_MINUTIAE; j++)
			{
				if(i != j && best_for_a_leg[i].idx_in_b == best_for_a_leg[j].idx_in_b 
					&& best_for_a_leg[i].idx_in_b == best_for_a[j].idx_in_b
					&& best_for_a_leg[j].supports > 0)
					flag = 1;
			}
			if(flag) continue;
			NumOfMatchedMinutiae++;
			MatchedListA[NumOfMatchedMinutiae] = i;
			MatchedListB[NumOfMatchedMinutiae * max_candidates] = best_for_a[i].idx_in_b;
			ScoreList[NumOfMatchedMinutiae * max_candidates] = corr_map[i][best_for_a[i].idx_in_b].accScore/corr_map[i][best_for_a[i].idx_in_b].hits;
			
			score += corr_map[i][best_for_a[i].idx_in_b].accScore;
			total_hits += corr_map[i][best_for_a[i].idx_in_b].hits;
			score += (corr_map_leg[i][best_for_a_leg[i].idx_in_b].accScore/2.0);
			total_hits += (corr_map_leg[i][best_for_a_leg[i].idx_in_b].hits/2);
		}
		*/
		// jea043005 --end--
	}
	NumOfMatchedMinutiae++;
#else /* defined jea050405 */
	NumOfMatchedMinutiae = 0;
	for(i = 0; i < CUBS_MAX_MINUTIAE; i++)
	{
		int flag = 0, k = 0;
		for(j = 0, k = j+1; j < CUBS_MAX_MINUTIAE && k < CUBS_MAX_MINUTIAE; k++)
		{
			if(corr_map[i][j].hits > 0)
			{
				flag = 1;
				if(corr_map[i][k].hits > 0
					&& corr_map[i][k].accScore/corr_map[i][k].hits < corr_map[i][j].accScore/corr_map[i][j].hits)
				{
					j = k;  // use the better one
				}
			}else
				j++;
		} // j
		if(flag)
		{
			MatchedListA[NumOfMatchedMinutiae] = i;
			MatchedListB[NumOfMatchedMinutiae * CUBS_MAX_CANDIDATES] = j;
			ScoreList[NumOfMatchedMinutiae * CUBS_MAX_CANDIDATES] = corr_map[i][j].accScore/corr_map[i][j].hits;
			NumOfMatchedMinutiae++;
		}
	}// i
#endif /* not define jea050405 */
		
	logfile_print("GetNumOfMatchedMinutiae(): ffter postprocess... got NumOfMatchedMinutiae=%d\n", NumOfMatchedMinutiae);

#ifdef CUBS_DEBUG
	fprintf(stderr,"After postprocess...\n");
	for(i = 0; i < NumOfMatchedMinutiae; i++){
		fprintf(stderr, "%d <-> %d (%8.6f)\n", MatchedListA[i], MatchedListB[i*max_candidates], ScoreList[i*max_candidates]);
	}
	fprintf(stderr,"Extend Match...\n");
#endif /* CUBS_DEBUG */
	// has problem if we use multiple choices here!!!!!!!!
	// added matched points to each others' neighbor list
#define jea050405
#ifndef jea050405
	for(i = 0; i < NumOfMatchedMinutiae; i++)
	{
		int idxCenterA = MatchedListA[i];
		int idxCenterB = MatchedListB[i*max_candidates];
		for(j = 0; j < NumOfMatchedMinutiae; j++)
		{
			if(j == i) continue;
			int k;
			PolarPt pt;
			pt.RefNo = MatchedListA[j];
			for(k=0; k < NeighborListA[idxCenterA].nCnt; k++)
			{
				if(pt.RefNo == NeighborListA[idxCenterA].arrNeighbors[k].RefNo)
					break;
			}
			if(k >= NeighborListA[idxCenterA].nCnt)
			{
				cart2polar(NeighborListA[pt.RefNo].m.m_nX, NeighborListA[pt.RefNo].m.m_nY, NeighborListA[idxCenterA].m.m_nX, NeighborListA[idxCenterA].m.m_nY, &(pt.radius), &(pt.angle));
				pt.angle = pt.angle - NeighborListA[idxCenterA].m.m_nTheta;
				if(pt.angle < 0) pt.angle += 360.0;
				pt.orientation = NeighborListA[pt.RefNo].m.m_nTheta - NeighborListA[idxCenterA].m.m_nTheta;
				if(pt.orientation < 0) pt.orientation += 360.0;
				(NeighborListA[idxCenterA].arrNeighbors).push_back(pt);
				NeighborListA[idxCenterA].nCnt++;
			}
			pt.RefNo = MatchedListB[j*max_candidates];
			for(k=0; k < NeighborListB[idxCenterB].nCnt; k++)
			{
				if(pt.RefNo == NeighborListB[idxCenterB].arrNeighbors[k].RefNo)
					break;
			}
			if(k >= NeighborListB[idxCenterB].nCnt)
			{
				cart2polar(NeighborListB[pt.RefNo].m.m_nX, NeighborListB[pt.RefNo].m.m_nY, NeighborListB[idxCenterB].m.m_nX, NeighborListB[idxCenterB].m.m_nY, &(pt.radius), &(pt.angle));
				pt.angle = pt.angle - NeighborListB[idxCenterB].m.m_nTheta;
				if(pt.angle < 0) pt.angle += 360.0;
				pt.orientation = NeighborListB[pt.RefNo].m.m_nTheta - NeighborListB[idxCenterB].m.m_nTheta;
				if(pt.orientation < 0) pt.orientation += 360.0;
				(NeighborListB[idxCenterB].arrNeighbors).push_back(pt);
				NeighborListB[idxCenterB].nCnt++;
			}
		} // j
	}//i
#else /* if define jea050405 */
	for(i = 0; i < NumOfMatchedMinutiae; i++)
	{
		int k = 0;
		PolarPt pt1, pt2;
		int idxCenterA = MatchedListA[i];
		for(j = 0; j < NumOfMatchedMinutiae; j++)
		{
			if(j == i) continue;
			pt1.RefNo = MatchedListA[j];
			for(k=0; k < NeighborListA[idxCenterA].nCnt; k++)
			{
				if(pt1.RefNo == NeighborListA[idxCenterA].arrNeighbors[k].RefNo)
					break;
			}
			if(k >= NeighborListA[idxCenterA].nCnt)
			{
				cart2polar(NeighborListA[pt1.RefNo].m.m_nX, NeighborListA[pt1.RefNo].m.m_nY, NeighborListA[idxCenterA].m.m_nX, NeighborListA[idxCenterA].m.m_nY, &(pt1.radius), &(pt1.angle));
				pt1.angle = pt1.angle - NeighborListA[idxCenterA].m.m_nTheta;
				if(pt1.angle < 0) pt1.angle += 360.0;
				pt1.orientation = NeighborListA[pt1.RefNo].m.m_nTheta - NeighborListA[idxCenterA].m.m_nTheta;
				if(pt1.orientation < 0) pt1.orientation += 360.0;
				(NeighborListA[idxCenterA].arrNeighbors).push_back(pt1);
				NeighborListA[idxCenterA].nCnt++;
			}
			for(int cnd1 = 0; cnd1 < CUBS_MAX_CANDIDATES; cnd1++)
			{
				int idxCenterB = MatchedListB[i*CUBS_MAX_CANDIDATES+cnd1];
				if(idxCenterB < 0)break;
				for(int cnd2 = 0; cnd2 < CUBS_MAX_CANDIDATES; cnd2++)
				{
					pt2.RefNo = MatchedListB[j*max_candidates+cnd2];
					if(pt2.RefNo < 0) break;
					for(k = 0; k < NeighborListB[idxCenterB].nCnt; k++)
					{
						if(pt2.RefNo == NeighborListB[idxCenterB].arrNeighbors[k].RefNo)
							break;
					}
					if(k >= NeighborListB[idxCenterB].nCnt)
					{
						cart2polar(NeighborListB[pt2.RefNo].m.m_nX, NeighborListB[pt2.RefNo].m.m_nY, NeighborListB[idxCenterB].m.m_nX, NeighborListB[idxCenterB].m.m_nY, &(pt2.radius), &(pt2.angle));
						pt2.angle = pt2.angle - NeighborListB[idxCenterB].m.m_nTheta;
						if(pt2.angle < 0) pt2.angle += 360.0;
						pt2.orientation = NeighborListB[pt2.RefNo].m.m_nTheta - NeighborListB[idxCenterB].m.m_nTheta;
						if(pt2.orientation < 0) pt2.orientation += 360.0;
						(NeighborListB[idxCenterB].arrNeighbors).push_back(pt2);
						NeighborListB[idxCenterB].nCnt++;
					}
				}//cnd2
			}//cnd1
		} // j
	}//i
#endif /* jea050405 */

	

	// perform extended match here.
	NumOfMatchedMinutiae = ExtendMatch(NeighborListA, /*cnt_a,*/
                NeighborListB, /*cnt_b,*/
                max_features, max_candidates,
                NumOfMatchedMinutiae,
                MatchedListA,    // size = max_features
                MatchedListB,    // size = max_features * max_candidates
                ScoreList       // size = max_features * max_candidates
                );
	

	// fill the result into pResultEx
	best_score = oo;
	score = 0.0;
	for(i = 0; i < NumOfMatchedMinutiae; i++)
	{
		pResultEx->FeaCorr[i].Lidx = MatchedListA[i];
		pResultEx->FeaCorr[i].Ridx = MatchedListB[i * max_candidates];
		pResultEx->FeaCorr[i].score = ScoreList[i * max_candidates];
		score += pResultEx->FeaCorr[i].score;
		if(ScoreList[i * max_candidates] < best_score)
		{
			best_score = ScoreList[i * max_candidates];
			best_idx   = i;
		}
	}
	// move the best matched pair to the front
	int tmp_int;
	double tmp_double;
	if(best_idx >= 0)
	{
		tmp_int = MatchedListA[0];
		MatchedListA[0] = MatchedListA[best_idx];
		MatchedListA[best_idx] = tmp_int;
		tmp_int = MatchedListB[0];
		MatchedListB[0] = MatchedListB[best_idx * max_candidates];
		MatchedListB[best_idx * max_candidates] = tmp_int;
		tmp_int = pResultEx->FeaCorr[0].Ridx;
		pResultEx->FeaCorr[0].Ridx = pResultEx->FeaCorr[best_idx].Ridx;
		pResultEx->FeaCorr[best_idx].Ridx = tmp_int;
		tmp_int = pResultEx->FeaCorr[0].Lidx;
		pResultEx->FeaCorr[0].Lidx = pResultEx->FeaCorr[best_idx].Lidx;
		pResultEx->FeaCorr[best_idx].Lidx = tmp_int;
		tmp_double = pResultEx->FeaCorr[0].score;
		pResultEx->FeaCorr[0].score = pResultEx->FeaCorr[best_idx].score;
		pResultEx->FeaCorr[best_idx].score = tmp_double;
	}
	// convert score into interval 0 and 1
	if(NumOfMatchedMinutiae > 0)
		ScoreList[0] = (NumOfMatchedMinutiae * MatchingParams.dT_MaxSFVDist - score)/(NumOfMatchedMinutiae * MatchingParams.dT_MaxSFVDist);
		//ScoreList[0] = (total_hits * MatchingParams.dT_MaxSFVDist - score)/(total_hits * MatchingParams.dT_MaxSFVDist);
	else
		ScoreList[0] = 0.0;

	pResultEx->NumOfMatchedMinutiae = NumOfMatchedMinutiae;
	
//#ifdef SHOW_MATCH_PTS
	for(i = 0; i < NumOfMatchedMinutiae; i++){
		//fprintf(stderr, "%3d <-> %3d\n", MatchedListA[i], MatchedListB[i*max_candidates]);
	}
//#endif /* SHOW_MATCH_PTS */

   return NumOfMatchedMinutiae;
}
