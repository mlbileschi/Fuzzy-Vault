#ifndef __COMMON_H__
#define __COMMON_H__
#pragma warning(disable: 4786)
#include    <vector>
#include    <map>
#include    <stdlib.h>
#include    <string>

#ifndef TRUE
#define TRUE        0x01
#endif
#ifndef FALSE
#define FALSE       0x00
#endif

#define CUBS_MAX_MINUTIAE  512

//#ifndef MAX_IMAGE_SIZE	
//#define MAX_IMAGE_SIZE	640
//#endif

/*
    File  : Common.h 
*/

/*--------------------------------------------------------
*Error constants
*---------------------------------------------------------*/
#define ERR_FILE_ERROR              0xF00
#define ERR_MEMORY_ERROR            0xF01
#define ERR_INVALID_PROPERTY        0xF02
#define ERR_PROPERTY_ALREADY_EXISTS 0xF03
#define ERR_INVALID_TYPE			0xF05
#define ERR_NO_ERROR                0x000

#define ERR_MISSING_IMAGE_INFO      0xE00
#define ERR_MISSING_IMAGE_DATA      0xE01
#define ERR_MISSING_PARAMETERS      0xE02

#define ERR_MISSING_MINUTIAE        0xD00
#define ERR_MISSING_RESULT          0xD01
#define ERR_TOO_MANY_FEATURES       0xD02

/*!
	@enum		DataType
	@abstract	Defines the data type of the property
*/
/*
enum DataType{		TYPE_UNKNOWN,
					TYPE_UINT8,
					TYPE_INT,
					TYPE_DOUBLE,
					TYPE_CHAR,
					TYPE_BINARY,
					TYPE_IMAGE
};
*/
/*!
    @struct Image
    @discussion 
        Contains the image information. It is assumed that the image is uint8 unless
		m_type has TYPE_DOUBLE
*/
/*
typedef struct
{
	int			m_nHeight;      //height
    int			m_nWidth;       //width
	void *		m_pData;		//image bits
	DataType	m_type;		    //defines data type
}Image;
*/
/*!
    @struct Minutiae
    @discussion
        Contains the minutiae descriptions
*/
typedef struct
{
	int		m_nFeatureNo;   //index of this feature
    int     m_nX;           //x -coordinate
    int     m_nY;           //y -coordinate
    int     m_nTheta;       //orientation
    char    m_cType;        //minutiae type
    int     m_nQuality;     //minutiae quality
}Minutiae;  // Used in CUBSEnroll and CUBSMatch


struct FPTemplate
{
	Minutiae marr[CUBS_MAX_MINUTIAE];
	int		 nCnt;
};
typedef struct FPTemplate* ptrFPTemplate;


/*!
  @struct MatchResult
  @discussion
     Contains the result of matching
*/
typedef struct
{
  bool   MatchPerformed;
  double similarity;
}MatchResult;

/*!
  @struct FeaPair
  @discussion
	Contains the feature pair that matched, and their matching score
*/
typedef struct
{
	int		Lidx; // feature number in first list
	int		Ridx; // feature number in second list
	double	score; // matching score
}FeaPair;

/*!
  @struct MatchResultEx
  @discussion
     Contains the extended matching result, i.e. the corresponding between minutiae numbers
*/
/*
typedef struct
{
	int     NumOfMatchedMinutiae;
	FeaPair FeaCorr[CUBS_MAX_MINUTIAE]; // Feature correspondence 
}MatchResultEx;
*/
typedef struct
{
	int     NumOfMatchedMinutiae;
	double	DistScore;
	double	PointScore;
	int		NumOfOverlap_a; // number of overlapped minutiae on A
	int		NumOfOverlap_b; // number of overlapped minutiae on B
	int		NumOfMinutiae_a;  // number of minutiae on A
	int		NumOfMinutiae_b;  // number of minutiae on B
	int     Width_a;          // width of print A (input)
	int     Height_a;         // height of print A (input)
	int     Width_b;          // width of print B (reference)
	int     Height_b;         // height of print B (reference)
	int     Width_comb;       // width of combined print
	int     Height_comb;      // height of combined print
	FeaPair FeaCorr[CUBS_MAX_MINUTIAE]; // Feature correspondence 
}MatchResultEx;

typedef MatchResultEx*       HCUBSMatchResultEx;
					
/*!
    @function UTL_ProcessFileName
    @abstract Conversts // to \\ (Cygwin)
*/
//void UTL_ProcessFileName(const char* pcInput,char* pcOutput);


/*!
    These are used to implement virtual mirror padding
*/
#define LEFT_MIR(x,wt)  (abs(x))
#define RIGHT_MIR(x,wt) ((2*(wt)-2)-(x))
#define TOP_MIR(y,ht)   (abs(y))
#define BOT_MIR(y,ht)   ((2*(ht)-2)-(y))

/*!
    @function PadWidth
    @discussion
        Implements mirror padding widthwise
*/
inline int PadWidth(int x,int wt)
{
    if(x < wt && x>=0)
        return x;
    if(x < 0)
        return LEFT_MIR(x,wt);
    return RIGHT_MIR(x,wt);
}

/*!
    @function PadHeight
    @discussion
        Implements mirror padding height wise
*/
inline int PadHeight(int y,int ht)
{
    if(y < ht && y>=0)
        return y;
    if(y <0)
        return TOP_MIR(y,ht);
    return BOT_MIR(y,ht);
}

#ifndef PI
#define PI              (3.1415926)
#endif

//Image* UTL_MakeImage(void *pData,int nHeight,int nWidth,DataType eType=TYPE_UINT8);

#endif //__COMMON_H__
