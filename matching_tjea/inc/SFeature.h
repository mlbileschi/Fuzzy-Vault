/*
  File : SFeature.h
*/
#ifndef __SFEATURE_H__
#define __SFEATURE_H__
#include <vector>
//#include    <map>
#include "../inc/common.h"

//#include <utility>
//using namespace std;

#define  ABSI(x) ((x)>0)?(x):(-1*(x))
#define  ABSD(x) ( ((x)>0.0) ? (x) : ((-1.0)*(x)) )
//#define  MAX(x,y) ((x>y)?x:y)

typedef unsigned char  QUADRANT; 
#define QUADRANT_0 0
#define QUADRANT_1 1
#define QUADRANT_2 2
#define QUADRANT_3 3
#define QUADRANT_4 4
#define QUADRANT_5 5
#define QUADRANT_6 6
#define QUADRANT_7 7

typedef std::pair<QUADRANT, QUADRANT> QPAIR;  // Quadrant pair for Secondary feature
typedef std::vector<int> SFVIDXVEC;           // Secondary feature pointers vector
typedef std::map<QPAIR, SFVIDXVEC> QUADSFVMAP;// Quadrant secondary features vector

/*!
  @struct PolarPt
  @discussion
    Contains the information of a point in terms of polar coordinates
*/
typedef struct {
	int    RefNo;            // Reference number refers back to minutia number.
	double radius;           // radius of the polar coordinates
	double angle;            // angle of the polar coordinates
	int    orientation;      // relative orientation w.r.t the reference pt
}PolarPt;

/*!
  @struct CartPt
  @discussion
    Contains the information of a point in terms of Cartesian coordinates
*/
typedef struct {
    int idx;
    int x;  // x coordiante
    int y;  // y coordinate
}CartPt;

/*!
  @struct SFeature
  @discussion
    Contains the information of a feature point in terms of secondary feature
*/
typedef struct {
  int          FeatureNo;   // index of current feature
  int          RefNo;       // Reference number refers back to minutia number.
  int          x;           // x coordiante
  int          y;           // y coordiante
  int          orientation; // orientation
  PolarPt      A;           // polar coordinates of the 1st neighbor w.r.t current pt
  PolarPt      B;           // polar coordinates of the 2nd neighbor w.r.t current pt
  double       angle;       // angle between the two legs
  int          ra;          // Relative orientation for 1st neighbor
  int          rb;          // Relative orientation for 2nd neighbor
  char         switched;    // Flag
#ifdef FV_IDX
  QPAIR        Q_Idx;       // Quadrant index (pair)
#endif
} SFeature, *SFeaturePtr;

// for internal using
typedef std::vector<PolarPt> PtVector;


typedef struct {
	Minutiae m;    // minutia 
	int		 nCnt; // number of neighbors
	PtVector arrNeighbors; // associate array of the neighbors
} NeighborInfo;

/*!
  @function cart2polar
  @abstract Convert a point from Cartesian coordinates into polar coordinates
  @param    cx    x coordinate of current point
  @param    cy    y coordiante of current point
  @param    rx    x coordinate of reference point
  @param    ry    y coordiante of reference point
  @param    r     radius of polar coordiante
  @param    theta angle of polar coordiante
*/
void cart2polar(int cx, int cy, int rx, int ry, double* r, double* theta);

/*!
  @function Minutiae2SFeature
  @abstract Convert Minutiae array into Secondary Feature Vectors array
  @param    nCnt       Number of features (minutiae/SFeature)
  @param    ArrM       Array of Minutiae
  @param    arrSFV     Array of SFeature
*/
#ifdef FV_IDX
	int Minutiae2SFeature(int nMCnt, Minutiae* arrM, int N, int mode, int* pnSFVCnt, SFeature** ppSFV, NeighborInfo* pNeighborList, QUADSFVMAP& QuadSFVList);
#else
	int Minutiae2SFeature(int nMCnt, Minutiae* arrM, int N, int mode, int* pnSFVCnt, SFeature** ppSFV, NeighborInfo* pNeighborList);
#endif
/*!
  @function ConvertFeaturePt
  @abstract Convert reference and test SFeature arrays into polar point arrays w.r.t the best-fit-points
  @param    list_a        Reference SFeature vectors array
  @param    cnt_a         Number of elements in list_a
  @param    best_idx_a    Index of best-fit-point in list_a
  @param    list_b        Test SFeature vectors array
  @param    cnt_b         Number of elements in list_b
  @param    best_idx_b    Index of best-fit-point in list_b
  @param    polar_list_a  Polar-points array of list_a
  @param    polar_list_b  Polar-points array of list_b
  @param    primary_od    Primary orientation difference between Reference and Test templates
*/
int ConvertFeaturePt(const SFeature *list_a, const int cnt_a, const int best_idx_a, const SFeature *list_b, const int cnt_b, const int best_idx_b, PolarPt* polar_list_a, PolarPt* polar_list_b, double primary_od);

#ifndef FN
/*!
  @function PolarListMatch
  @abstract Find the one-to-one correspondence between two polar-point lists
  @param    polar_list_a  The 1st polar-point list
  @param    cnt_a         Number of elements in polar_list_a
  @param    polar_list_b  The 2nd polar-point list
  @param    cnt_b         Number of elements in polar_list_b
  @param    MatchedListA  Matched point indices array of polar_list_a (at least has the same size as polar_list_a)
  @param    MatchedListB  Matched point indices array of polar_list_b (at least has the same size as polar_list_b)
  @param    score         Total distance of matched pairs.
  @param    mode          Regular or Forced mode
*/
int PolarListMatch(const PolarPt* polar_list_a, const int cnt_a, const PolarPt* polar_list_b, const int cnt_b, int *MatchedListA, int *MatchedListB, double* score, int mode);
#else /* def FN */
/*
** Function: PolarFNMatch
** Points matching by using minimum cost flow.
** Returns:  the number of matched polar pts.
**           -1 if anything wrong.
** Parameters:
**   [IN] polar_list_a  // 1st polar points list
**   [IN] cnt_a         // number of points in 1st polar points list
**   [IN] polar_list_b  // 2nd polar points list
**   [IN] cnt_b         // number of points in 2nd polar points list
**   [OUT] MatchedListA // should be allocated outside
**   [OUT] MatchedListB // should be allocated outside
**   [out] score        // final total score
*/
int
PolarFNMatch(
			 PolarPt* polar_list_a, 
			 const int cnt_a, 
			 PolarPt* polar_list_b, 
			 const int cnt_b,
			 int *MatchedListA, // same size as polar_list_a
			 int *MatchedListB, // same size as polar_list_b
			 double *score
			 );
#endif /* FN */
/*!
  @function FeatureMatch
  @abstract Match two SFeature vectors array
  @param    list_a          1st SFeature vectors array
  @param    list_a_sz       Number of elements in list_a
  @param    list_b          2nd SFeature vectors array
  @param    list_b_sz       Number of elements in list_b
  @param    max_features    Maximum number of features
  @param    max_candidates  Maximum number of matched candidates for each feature
  @param    NumMatched      Number of matched feature points
  @param    MatchedListA    Array of indices of matched features in list_a
  @param    MatchedListB    Array of indices of matched features in list_b
  @param    ScoreList       Array of scores for every matched feature pairs
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
             );

/*!
  @function PostProcess
  @abstract Filterint out false matched SFeature pairs and estimating primary orientation difference between two SFeature vectors lists
  @param    list_a          1st SFeature vectors array
  @param    list_b          2nd SFeature vectors array
  @param    max_features    Maximum number of features
  @param    max_candidates  Maximum number of matched candidates for each feature
  @param    NumMatched      Number of matched SFeature pairs
  @param    MatchedListA    Array of indices of matched features in list_a
  @param    MatchedListB    Array of indices of matched features in list_b
  @param    ScoreList       Array of scores for every matched feature pairs
*/
double PostProcess(
                   const SFeature *list_a,
				   const SFeature *list_b,
				   const int      max_features,
                   const int      max_candidates,
                   int*           NumMatched,
                   int*           MatchedListA,    // size = max_features
                   int*           MatchedListB,    // size = max_features * max_candidates
                   double         *ScoreList       // size = max_features * max_candidates
                   );

/*!
  @function GetNumOfMatchedMinutiae
  @abstract Getting number of matched minutiae from two SFeature vectors arrays
  @param    list_a          1st SFeature vectors array
  @param    list_a_sz       Number of elements in list_a
  @param    list_b          2nd SFeature vectors array
  @param    list_b_sz       Number of elements in list_b
  @param    max_features    Maximum number of features
  @param    max_candidates  Maximum number of matched candidates for each feature
  @param    NumMatched      Number of matched feature points
  @param    MatchedListA    Array of indices of matched features in list_a
  @param    MatchedListB    Array of indices of matched features in list_b
  @param    ScoreList       Array of scores for every matched feature pairs
  @param    primary_od      Primary orientation difference between list_a and list_b
  @param    mode            Specify the matching mode (regular or forced)
*/
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
                            );

/*!
  @function SFV_VerifyMatch
  @abstract Top-level functino to match two SFeature vectors arrays
  @param    sfv_arr           1st (Test) SFeature vectors array
  @param    sfv_cnt           Number of elements in sfv_arr
  @param    Temp_sfv_arr      2nd (Reference) SFeature vectors array
  @param    Temp_sfv_cnt      Number of elements in Temp_sfv_arr
  @param    MatchingPerformed Indicates if matching is performed
  @param    similarity        Similarity level of Test and Reference arrays (between 0 and 1)
  @param    pResultEx         Contains the pairing information of matched feature points
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
  std::vector<std::pair<int,int> > &matchingPairs);
/*
int SFV_VerifyMatch(
					int nMCntTest, Minutiae* arrMTest,
					int nMCntRef, Minutiae* arrMRef,
                    SFeature         *sfv_arr, 
                    const int        sfv_cnt,
                    const SFeature   *Temp_sfv_arr, 
                    const int        Temp_sfv_cnt,
                    bool             *MatchingPerformed,
                    double           *similarity,
					MatchResultEx    *pResultEx
                    );
*/
/*!
  @function NumMinutiaeInCH
  @abstract Calaulate the number of minutiae in the Convex Hull of Matched minutiae set
  @param    sfv     SFeature vectors array
  @param    SfvCnt  Number of elements in sfv
  @param    MatchedList Array of indeices of matched minuatiae points
  @param    MatchedCnt  Number of matched points
  @param    BB_minx     Minimum x value of the bounding box of the Hull
  @param    BB_maxx     Maximum x value of the bounding box of the Hull
  @param    BB_miny     Minimum y value of the bounding box of the Hull
  @param    BB_maxy     Maximum y value of the bounding box of the Hull
 */
int NumMinutiaeInCH(const SFeature *sfv, const int SfvCnt, const int* MatchedList, const int MatchedCnt,
                    int* BB_minx, int* BB_miny, int* BB_maxx, int* BB_maxy);

/*!
  @function NumMinutiaeInCH2
  @abstract Calculate the number of minutiae of print A that are within the convex hull of print B
  @param    sfv_A     SFeature vectors array of print A
  @param    SfvCnt_A  Number of elements in sfv_A
  @param    sfv_B     SFeature vectors array of print B
  @param    SfvCnt_B  Number of elements in sfv_B
*/
int NumMinutiaeInCH2(
                     const SFeature *sfv_A, const int SfvCnt_A, 
                     const SFeature *sfv_B, const int SfvCnt_B
                     );

int NumMinutiaeInCH3(const Minutiae *arrM_A, const int Cnt_A, 
                     const Minutiae *arrM_B, const int Cnt_B,
					 int* BB_minx, int* BB_miny, int* BB_maxx, int* BB_maxy
                     );

/*! 
  @function NumOverlappedMinutiae
*/
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
                          );
#endif /* noNEWSFV */
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
						  MatchResultEx    *pResultEx
                          );

int ExtendMatch(NeighborInfo* NeighborListA,/* const int cnt_a,*/
                NeighborInfo* NeighborListB,/* const int cnt_b,*/
                const int max_features,
                const int max_candidates,
                const int NumMatched,
                int*      MatchedListA,    // size = max_features
                int*      MatchedListB,    // size = max_features * max_candidates
                double    *ScoreList       // size = max_features * max_candidates
				);
/// Utility functions ///
double abs_angle_diff(double x, double y); // defined in polar_util.cpp
#endif //__SFEATURE_H__
