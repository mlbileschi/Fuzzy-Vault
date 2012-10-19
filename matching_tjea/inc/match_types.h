#ifndef __MATCH_TYPES_H__
#define __MATCH_TYPES_H__

#include "../inc/common.h"

/*
	algorithm parameters
*/
#define MAX_MINUTIAE	200
#define MAX_PAIRS		100
#define MAX_TPS_ORDER	15
#define TPS_PAIRS		6
/*
	error codes
*/
#define UTL_ERR_NO_ERROR			0
#define UTL_ERR_INVALID_PARAMETER	1
#define UTL_ERR_INTERNAL_ERROR		2
#define UTL_ERR_INSUFFICIENT_POINTS 3
/*
	radial match parameters
*/
#define DR			(0.05)	
#define DTHETA		(10)
#define DDELTA		(20)
/*
final match parameters
*/
#define  YTOL	6
#define	 XTOL	6
#define  RTOL	8

#if defined(__cplusplus)
extern "C" {
#endif


typedef Minutiae Minutia;

/*!
	@struct RMinutia
*/
typedef struct __tagRMinutia
{
	double	m_dR;
	int		m_nTheta;
	int		m_nDelta;
}RMinutia;

/*!
	@struct Pair
*/
typedef struct __tagPair
{
	int		m_nRIdx;  /*reference minutia index*/
	int		m_nTIdx;  /*test minutia index*/
	double  m_dScore; /*pair match score*/
}Pair;

/*!
    @struct   TxParams
	@abstract contains Thin plate spline parameters
*/
typedef struct __tagTxParams
{
	int		m_nCnt;	
	double  m_dTxX[MAX_PAIRS];
	double  m_dTxY[MAX_PAIRS];
	Minutia m_minCtrl[MAX_MINUTIAE]; //control points
}TxParams;

#if defined(__cplusplus)
}
#endif

#endif
