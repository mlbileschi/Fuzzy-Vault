#include "../inc/global_parameters.h"

std::vector<double> Global_Parameters::DISTRO = std::vector<double>(150);

// Algorithm choice parameters
FEATURE_EXTRACTION_ALGORITHM Global_Parameters::FEATURE_EXTRACTION;	
ENHANCEMENT_ALGORITHM Global_Parameters::ENHANCEMENT;	
MATCHING_ALGORITHM Global_Parameters::MATCHING;	
// Display and debug parameters
bool Global_Parameters::SAVE_IMAGES;
bool Global_Parameters::PRINT_LOGFILE;
//----------------------------
//	Enhancement
//----------------------------
const double Global_Parameters::E_ALPHA =	0.5;	//Pseudo filtering
const double Global_Parameters::E_FLOW =			5;	//smallest ridge distance
const double Global_Parameters::E_FHIGH =		15;	//largest ridge distance
const int Global_Parameters::E_BLKSZ=			12;	//block size
const int Global_Parameters::E_OVRLP=			6;	//size of overlap
const double Global_Parameters::E_STRETCH	=	10;	//contrast stretch
const double Global_Parameters::E_THRESH	=	10;	//threshold for energy map  ////// CHECK THIS!!!!!
const double Global_Parameters::E_NEGLECT_ENERGY_MAP=	1;	//don't use energy map
const double Global_Parameters::E_PREFILTERING	=	1;	//use pseudo matched filter
const double Global_Parameters::E_CONTEXTFILTERING	=1;	//do contextual filtering
const double Global_Parameters::E_FFTN	=32;	
const double Global_Parameters::E_FILTORD	=8;	
const double Global_Parameters::E_FRAC	=0.89;	
const double Global_Parameters::E_STRICT =1;	
//const unsigned char  Global_Parameters::E_FILL =254; //	#average background pixel (should be changed for different databases)
int Global_Parameters::E_PAD_XL;	
int Global_Parameters::E_PAD_XR;	
int Global_Parameters::E_PAD_YL;	
int Global_Parameters::E_PAD_YR;	
//----------------------------
//	Feature Extraction
//----------------------------
const double Global_Parameters::F_BLKSZ  =       	12;
const double Global_Parameters::F_THRESH    =    	8;
const double Global_Parameters::F_ENERGYTH  =    	27;
const double Global_Parameters::F_SEG       =    	4;	
//#----------------------------
//#	Matching parameters
//#----------------------------
const double Global_Parameters::M_W_RADIUS	=	500.0;
const double Global_Parameters::M_W_ORIENTATION=		1.0;
const double Global_Parameters::M_W_ANGLE	=	1.0;
const double Global_Parameters::M_T_ORIENTATION=	  	30;
const double Global_Parameters::M_T_P_ORIENTATION =	20;
const double Global_Parameters::M_T_RADIUS	=	0.15;
const double Global_Parameters::M_T_P_RADIUS	=	0.1;
const double Global_Parameters::M_T_MINANGLETOL=		15;
const double Global_Parameters::M_T_P_MINANGLETOL=	10;
const double Global_Parameters::M_T_MAXSFVDIST	=	4.5;
const double Global_Parameters::M_T_MINSFVMATCH=		1;
const double Global_Parameters::M_F_HYBRID	=	1;
const double Global_Parameters::M_T_F_ORIENTATION =	20;
const double Global_Parameters::M_T_F_RADIUS=		0.10;
const double Global_Parameters::M_T_F_MINANGLETOL	=5;
const double Global_Parameters::M_T_MINSFVREQ       =    25;
const double Global_Parameters::M_T_TRIES         =      1;
//#---------------------------
//# New matching parameters
//#---------------------------
const double Global_Parameters::M_T_TRYNREF     =        10 ; //try first n matched sfv as reference points
const double Global_Parameters::M_T_GOODWIDTH   =        400;
const double Global_Parameters::M_T_GOODHEIGHT   =       400;
const double Global_Parameters::M_T_F_MAXRADIUS  =       2.0;
const double Global_Parameters::M_T_F_MAXANGLE     =     2.0;
const double Global_Parameters::M_T_F_MAXORIENTATION  =  2.0;
const double Global_Parameters::M_T_P_MAXRADIUS      =   1.0;
const double Global_Parameters::M_T_P_MAXANGLE       =   2.0;
const double Global_Parameters::M_T_P_MAXORIENTATION  =  2.0;
const double Global_Parameters::M_F_BRUTEFORCE        =  0;

	
void Global_Parameters::SET_CURRENT_CONFIGURATION(DATABASE_CONFIGURATION)
{

	return;
}
	
	/*
class Global_Parameters{
public:
	//----------------------------
	//	Enhancement
	//----------------------------
	static const double E_ALPHA =	0.5;	//Pseudo filtering
	static const double E_FLOW =			5;	//smallest ridge distance
	static const double E_FHIGH =		15;	//largest ridge distance
	static const double E_BLKSZ=			12;	//block size
	static const double E_OVRLP=			6;	//size of overlap
	static const double E_STRETCH	=	10;	//contrast stretch
	static const double E_THRESH	=	27;	//threshold for energy map
	static const double E_NEGLECT_ENERGY_MAP=	1;	//don't use energy map
	static const double E_PREFILTERING	=	1;	//use pseudo matched filter
	static const double E_CONTEXTFILTERING	=1;	//do contextual filtering
	//----------------------------
	//	Feature Extraction
	//----------------------------
	static const double F_BLKSZ  =       	12;
	static const double F_THRESH    =    	8;
	static const double F_ENERGYTH  =    	27;
	static const double F_SEG       =    	4;	
	//#----------------------------
	//#	Matching parameters
	//#----------------------------
	static const double M_W_RADIUS	=	500.0;
	static const double M_W_ORIENTATION=		1.0;
	static const double M_W_ANGLE	=	1.0;
	static const double M_T_ORIENTATION=	  	30;
	static const double M_T_P_ORIENTATION =	20;
	static const double M_T_RADIUS	=	0.15;
	static const double M_T_P_RADIUS	=	0.1;
	static const double M_T_MINANGLETOL=		15;
	static const double M_T_P_MINANGLETOL=	10;
	static const double M_T_MAXSFVDIST	=	4.5;
	static const double M_T_MINSFVMATCH=		1;
	static const double M_F_HYBRID	=	1;
	static const double M_T_F_ORIENTATION =	20;
	static const double M_T_F_RADIUS=		0.10;
	static const double M_T_F_MINANGLETOL	=5;
	static const double M_T_MINSFVREQ       =    25;
	static const double M_T_TRIES         =      1;
	//#---------------------------
	//# New matching parameters
	//#---------------------------
	static const double M_T_TRYNREF     =        10 ; //try first n matched sfv as reference points
	static const double M_T_GOODWIDTH   =        400;
	static const double M_T_GOODHEIGHT   =       400;
	static const double M_T_F_MAXRADIUS  =       2.0;
	static const double M_T_F_MAXANGLE     =     2.0;
	static const double M_T_F_MAXORIENTATION  =  2.0;
	static const double M_T_P_MAXRADIUS      =   1.0;
	static const double M_T_P_MAXANGLE       =   2.0;
	static const double M_T_P_MAXORIENTATION  =  2.0;
	static const double M_F_BRUTEFORCE        =  0;
};
*/
