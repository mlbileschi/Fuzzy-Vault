//#include "stdafx.h"
#include <stdio.h>
#include <stdarg.h>
#include "../inc/Global_Parameters.h"

void logfile_print(char* pFormat, ...)
{
	if(!Global_Parameters::PRINT_LOGFILE)
		return;
//	CTime t = CTime::GetCurrentTime();
	
//	CString s = t.Format("[%Y/%m/%d %H:%M:%S]");
	char    chMsg[4096];
	char    tmpMsg[4096];
	
//	_tcsncpy(chMsg, (LPCTSTR)s, s.GetLength());
	
	va_list pArg;
	
	va_start(pArg, pFormat);
	vsprintf(tmpMsg, pFormat, pArg);
	va_end(pArg);
	
	::sprintf(chMsg,"%s", tmpMsg);
	
	//_putts(chMsg);
	FILE *logfile=fopen("CUBS.log", "a+");
	if(logfile!=NULL) {
		fprintf(logfile,chMsg);
		fclose(logfile);	
	}
}

