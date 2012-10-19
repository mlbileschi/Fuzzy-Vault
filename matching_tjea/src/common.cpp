/*!
    @file    common.cpp
    @discussion
            Provides implementation of the common functions
*/
#include "../inc/common.h"

//-------------------------------------------------------------------
//UTL_ProcessFileName
//-------------------------------------------------------------------
void UTL_ProcessFileName(const char* pcInput,char* pcOutput)
{
	while(*pcOutput=*pcInput)
	{
		if(*pcOutput=='/')
			*pcOutput='\\';
		pcInput++;
		pcOutput++;
	}
}



//-------------------------------------------------------------------
//UTL_MakeImage
//-------------------------------------------------------------------
/*
Image* UTL_MakeImage(void * pData,int nHeight,int nWidth,DataType eType)
{
	Image	*img	=	(Image*)malloc(sizeof(Image)); //XIA
	if(img == NULL)
		return NULL;
	img->m_pData	=	pData;
	img->m_nWidth	=	nWidth;
	img->m_nHeight	=	nHeight;
	img->m_type		=	eType;
	return img;
}
*/
