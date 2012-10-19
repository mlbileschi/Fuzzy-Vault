#include <stdio.h>
#include <stdlib.h>
#include "../inc/SFeature.h"
#include "../inc/matching.h"
#include <vector>
#include <queue>
#include <map>

using namespace std;
typedef queue<int> IDXQUEUE;
typedef queue<FeaPair> MATCHQUEUE;

// copied from polar_util.cpp
static double radius_dist(const PolarPt pt1, const PolarPt pt2)
{
	double tolerate_ratio = 0.1;
	double unit = pt2.radius * tolerate_ratio;
	double dist;
	if(unit < 7.0) unit = 7.0;
	if(unit > 15.0) unit = 15.0; // jea050505
	//if(unit < 5.0) unit = 5.0;// jea043005
	if(unit)
		dist = (ABSD(pt1.radius - pt2.radius))/unit;
	else
		dist = 0.0;

	if(dist > MatchingParams.dT_P_MaxRadius) dist = 1000.0;

	return dist;
}

static double angular_dist(const PolarPt pt1, const PolarPt pt2)
{
	double max_r = 100.0;
	double min_r = 5.0;
	//jea050505 start
	double lb = 5.0;  // lower bound of tolerable angle difference when radius longer than max_r pixels
	double ub = 15.0; // upper bound of tolerable angle difference when radius shorter than min_r pixels //jea042905
	//double lb = 10.0;  // lower bound of tolerable angle difference when radius longer than max_r pixels
//	double ub = 20.0; // upper bound of tolerable angle difference when radius shorter than min_r pixels //jea042905
	//jea050505 end
	double dist = 0.0;
	double unit = 0.0;

	if(pt2.radius >= max_r)
		unit = lb;
	else if(pt2.radius <= min_r)
		unit = ub;
	else
		unit = ub - ((pt2.radius - min_r) * (ub - lb) / (max_r - min_r));

	dist = abs_angle_diff(pt1.angle, pt2.angle) / unit;

    if(dist > MatchingParams.dT_P_MaxAngle) dist = 1000.0;

	return dist;
}

static double orientation_dist(const PolarPt pt1, const PolarPt pt2)
{
	double max_r = 100.0;
	double min_r = 5.0;
	double lb = 10.0;  // lower bound of tolerable angle difference when radius longer than max_r pixels
	double ub = 30.0; // upper bound of tolerable angle difference when radius shorter than min_r pixels
	double dist = 0.0;
	double unit = 0.0;

	if(pt2.radius >= max_r)
		unit = lb;
	else if(pt2.radius <= min_r)
		unit = ub;
	else
		unit = ub - ((pt2.radius - min_r) * (ub - lb) / (max_r - min_r));

	dist = abs_angle_diff(pt1.orientation, pt2.orientation) / unit;

	if(dist > MatchingParams.dT_P_MaxOrientation) dist = 1000.0;

	return dist;
}

/*!
  @function: pure_polar_dist
  @discussion:
    Rreturn the distance between two features without thresholding
*/
static double pure_polar_dist(const PolarPt pt1, const PolarPt pt2)
{
	double dist = 0.0;
	//dist = radius_dist(pt1, pt2) + angular_dist(pt1, pt2) + orientation_dist(pt1, pt2);
	dist = 2.0 * radius_dist(pt1, pt2) + 0.5 * angular_dist(pt1, pt2) + 0.5*orientation_dist(pt1, pt2);
	
	return dist;
}
// end copied from polar_util.cpp

static int expend(IDXQUEUE& queueA, IDXQUEUE& queueB, NeighborInfo* NeighborListA, NeighborInfo* NeighborListB, MATCHQUEUE& queueMatch)
{
	bool marksA[CUBS_MAX_MINUTIAE]={false};
	bool marksB[CUBS_MAX_MINUTIAE]={false};
	int  nCnt = 0;
	int  i = 0, j = 0;

	if(queueA.empty() || queueB.empty())
		return nCnt;

	nCnt++;
	marksA[queueA.front()] = true;
	marksB[queueB.front()] = true;
	
	// match on neighbors
	while(!queueA.empty() && !queueB.empty())
	{
		FeaPair  pairMatch;
		int idxA = queueA.front();
		int idxB = queueB.front();
		queueA.pop();
		queueB.pop();
		for(i = 0; i < NeighborListA[idxA].nCnt; i++)
		{
			int     bestB = -1;
			double  best_dist = 10000.0;
			PolarPt ptA = NeighborListA[idxA].arrNeighbors[i];

			if(marksA[ptA.RefNo])
			{
				continue;
			}
			for( j = 0; j < NeighborListB[idxB].nCnt; j++)
			{
				double dist      = 10000.0;
				PolarPt ptB = NeighborListB[idxB].arrNeighbors[j];
				// if the angle difference in smaller than 30 degree
				if(marksB[ptB.RefNo])
				{
					continue;
				}
				if(ABSD(ptA.angle - ptB.angle) < 25.0 || ABSD(ptA.angle - ptB.angle) > 335.0)
				{
					dist = pure_polar_dist(ptA, ptB);
					if(dist < 3.5 && dist < best_dist)// the less the better
					//if(dist < 3.0 && dist < best_dist)// the less the better
					{
						best_dist = dist;
						bestB = ptB.RefNo;
					}
				}
			} // j
			if(bestB > -1)
			{
				pairMatch.Lidx  = ptA.RefNo;
				pairMatch.Ridx  = bestB;
				pairMatch.score = best_dist;
				// add into queue
				queueA.push(pairMatch.Lidx);
				queueB.push(pairMatch.Ridx);
				queueMatch.push(pairMatch);
				marksA[pairMatch.Lidx] = true;
				marksB[pairMatch.Ridx] = true;				
				int t = queueMatch.size();
				nCnt++;
			}
		} // i
	}
	return nCnt;
}

// return matched minutiae points
int ExtendMatch(NeighborInfo* NeighborListA,/* const int cnt_a,*/
                NeighborInfo* NeighborListB,/* const int cnt_b,*/
                const int max_features,
                const int max_candidates,
                const int NumMatched,
                int*      MatchedListA,    // size = max_features
                int*      MatchedListB,    // size = max_features * max_candidates
                double    *ScoreList       // size = max_features * max_candidates
				)
{
	MATCHQUEUE queueBestMatch;
	int      i = 0, j = 0;
	int		 NumExtendMatch = 0;
	int      CurrentBestNumExtendMatch = 0;
	double best_score = 10000; //jea042905
	map<int,int> mapMatch;
	int          cntCandidates = 0;
	
	for(i = 0; i < NumMatched; i++)
	{
		IDXQUEUE queueA;
		IDXQUEUE queueB;
		MATCHQUEUE queueMatch;
		FeaPair pairMatch;
		queueA.push(MatchedListA[i]);
		for(j = 0; j < CUBS_MAX_CANDIDATES; j++)
		{
			if(MatchedListB[i * max_candidates+j] < 0) break;
			queueB.push(MatchedListB[i * max_candidates+j]);
			pairMatch.Lidx = MatchedListA[i];
			pairMatch.Ridx = MatchedListB[i * max_candidates + j];
			pairMatch.score = ScoreList[i * max_candidates + j];
			queueMatch.push(pairMatch);
			if(mapMatch[pairMatch.Lidx] == pairMatch.Ridx)
			{
				continue;
			}
			else
			{
				cntCandidates++;
			}
			NumExtendMatch = expend(queueA, queueB, NeighborListA, NeighborListB, queueMatch);
			
			if(NumExtendMatch >= 3)
			{
				MATCHQUEUE tmpMQ;
				double scores = 0.0;
				tmpMQ = queueMatch;
				while(!tmpMQ.empty())
				{
					scores += tmpMQ.front().score;
					tmpMQ.pop();
				}
				scores = scores / (double)NumExtendMatch;
				
				if(NumExtendMatch > CurrentBestNumExtendMatch)
				{
					CurrentBestNumExtendMatch = NumExtendMatch;
					queueBestMatch = queueMatch;
					best_score = scores;
				}
				else if(NumExtendMatch == CurrentBestNumExtendMatch)
				{
					if(scores < best_score)
					{
						CurrentBestNumExtendMatch = NumExtendMatch;
						queueBestMatch = queueMatch;
						best_score = scores;
					}
				}
				if(NumExtendMatch >= CurrentBestNumExtendMatch)
				{
					mapMatch.clear();
					while(!queueMatch.empty())
					{
						int a = queueMatch.front().Lidx;
						int b = queueMatch.front().Ridx;
						mapMatch[a] = b;
						queueMatch.pop();
					}
				}
			}
		}// j
	} // i
	i = 0;
	unsigned char flag = 0; // possible imposter
	if(CurrentBestNumExtendMatch <= 6 && 2*(cntCandidates-1) > NumMatched)
	{
		flag = 1;
//#####		fprintf(stderr,"!!!!! possible imposter !!!!!\n");
	}
	while(!queueBestMatch.empty())
	{
		MatchedListA[i] = queueBestMatch.front().Lidx;
		MatchedListB[i * max_candidates] = queueBestMatch.front().Ridx;
		if(flag) // possible imposter
		{
			ScoreList[i * max_candidates] = queueBestMatch.front().score + 1.0;
		}
		else
		{
			ScoreList[i * max_candidates] = queueBestMatch.front().score;
		}
		queueBestMatch.pop();
		i++;
	}
	return CurrentBestNumExtendMatch;
}
