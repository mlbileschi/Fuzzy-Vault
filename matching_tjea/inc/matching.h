#ifndef __MATCHING_H__
#define __MATCHING_H__
//#include "SFeature.h"
#include <vector>
#include "../inc/common.h"

#define CUBS_MAX_FEATURES   5000//1024
#define CUBS_MAX_CANDIDATES 5
#define CUBS_NUM_OF_NEIGHBORS 6//6 jea043005
//#define TOLERANCE      0.10
#define NUM_OF_BINS    30 //jea050105
#define MODE_FORCED    1
#define MODE_REGULAR   2
#ifdef CHEATING
#define CHEATING_ANGLE 50.0
#endif

/* Error code for matching module */
#define MATCH_ERR_INVALID_INPUT  0xD10
#define MATCH_ERR_INTERNAL_ERROR 0xD11

/*!
  @struct Params
  @discussion
    Contains all the parameters/thresholds we need to perform the matching
*/
typedef struct {
  double dT_MaxSFVDist;		// Max. allowed distance between Secondary Feature Vectors
  double dT_MinSFVMatch;	// Min. required SVF match to claim a positive genuine match
  int    nF_Hybrid;         // Hybrid flag (combine regular SFV match and brute-force match)
  int    nT_MinSFVReq;      // Min. number of SFV required for non-brute-force matching (in hybrid mode)
  int    nT_GoodWidth;      // acceptable width after combining query and template images
  int    nT_GoodHeight;     // acceptable height after combining query and template images
  int    nT_TryNRef;        // To specify the number of matched SFVs that we use as ref. points
                            //  to find final number of matches. -1 means try all.
                            /* Following are the thresolds for new feature distance */
  double dT_F_MaxRadius;    // Max acceptable radius (in UNITS) for considering a SFV matching.
  double dT_P_MaxRadius;    // Max acceptable radius (in UNITS) for considering a point matching.
  double dT_F_MaxAngle;     // Max acceptable angle (in UNITS) for considering a SFV matching.
  double dT_P_MaxAngle;     // Max acceptable angle (in UNITS) for considering a point matching.
  double dT_F_MaxOrientation;  // Max acceptable Orientation (in UNITS) for considering a SFV matching.
  double dT_P_MaxOrientation;  // Max acceptable Orientation (in UNITS) for considering a point matching.
  int    nF_BruteForce;		// flag to specify if we are using brute force matching only
}Params;

extern Params MatchingParams;

int CUBS_MatchFeaturesInternal(ptrFPTemplate hRef, ptrFPTemplate hTst, 
							   MatchResult*	pMatchResult, MatchResultEx* pResultEx, std::vector<std::pair<int,int> > &matchingPairs);
int CUBS_Match_IBM_Hash(ptrFPTemplate hRef, ptrFPTemplate hTst, 
							   MatchResult*	pMatchResult, MatchResultEx* pResultEx);


#endif //__MATCHING_H__
