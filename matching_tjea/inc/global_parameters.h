#include <vector>


enum FEATURE_EXTRACTION_ALGORITHM {
	ANSI_NIST_MINDTCT_FEATURES, CUBS_FEATURES };
enum ENHANCEMENT_ALGORITHM {
	ENHANCEMENT_STFT, ENHANCEMENT_DMF, ENHANCEMENT_MRF };
enum MATCHING_ALGORITHM {
	MATCHING_MCF, MATCHING_GRAPH, MATCHING_SYMMETRIC_HASH, MATCHING_IBM_HASH, MATCHING_FUZZY_VAULT};
enum DATABASE_CONFIGURATION { CUSTOM,
	FVC2002_DB1, FVC2002_DB2, FVC2002_DB3, FVC2002_DB4, 
	FVC2004_DB1, FVC2004_DB2, FVC2004_DB3, FVC2004_DB4, 
	FVC2006_DB1, FVC2006_DB2, FVC2006_DB3, FVC2006_DB4};

class Global_Parameters{
public:
	static std::vector<double> DISTRO;
	// Algorithm choice parameters
	static FEATURE_EXTRACTION_ALGORITHM FEATURE_EXTRACTION;	
	static ENHANCEMENT_ALGORITHM ENHANCEMENT;	
	static MATCHING_ALGORITHM MATCHING;	
	// Adjust matching parameters depending on set database configuration
	static void SET_CURRENT_CONFIGURATION(DATABASE_CONFIGURATION);
	// Display and debug parameters
	static bool SAVE_IMAGES;
	static bool PRINT_LOGFILE;
	//----------------------------
	//	Enhancement
	//----------------------------
	static const double E_ALPHA ;	//Pseudo filtering
	static const double E_FLOW ;	//smallest ridge distance
	static const double E_FHIGH ;	//largest ridge distance
	static const int E_BLKSZ ;	//block size
	static const int E_OVRLP ;	//size of overlap
	static const double E_STRETCH;	//contrast stretch
	static const double E_THRESH;	//threshold for energy map
	static const double E_NEGLECT_ENERGY_MAP;	//don't use energy map
	static const double E_PREFILTERING;	//use pseudo matched filter
	static const double E_CONTEXTFILTERING;	//do contextual filtering
	static const double E_FFTN;	
	static const double E_FILTORD;	
	static const double E_FRAC;	
	static int E_PAD_XL;	
	static int E_PAD_XR;	
	static int E_PAD_YL;	
	static int E_PAD_YR;	
	static const double E_STRICT;	
//	static const unsigned char E_FILL;	
	//----------------------------
	//	Feature Extraction
	//----------------------------
	static const double F_BLKSZ ;
	static const double F_THRESH ;
	static const double F_ENERGYTH ;
	static const double F_SEG ;	
	//#----------------------------
	//#	Matching parameters
	//#----------------------------
	static const double M_W_RADIUS;
	static const double M_W_ORIENTATION;
	static const double M_W_ANGLE;
	static const double M_T_ORIENTATION;
	static const double M_T_P_ORIENTATION ;
	static const double M_T_RADIUS;
	static const double M_T_P_RADIUS;
	static const double M_T_MINANGLETOL;
	static const double M_T_P_MINANGLETOL;
	static const double M_T_MAXSFVDIST;
	static const double M_T_MINSFVMATCH;
	static const double M_F_HYBRID;
	static const double M_T_F_ORIENTATION ;
	static const double M_T_F_RADIUS;
	static const double M_T_F_MINANGLETOL;
	static const double M_T_MINSFVREQ ;
	static const double M_T_TRIES ;
	//#---------------------------
	//# New matching parameters
	//#---------------------------
	static const double M_T_TRYNREF  ; //try first n matched sfv as reference points
	static const double M_T_GOODWIDTH ;
	static const double M_T_GOODHEIGHT  ;
	static const double M_T_F_MAXRADIUS ;
	static const double M_T_F_MAXANGLE  ;
	static const double M_T_F_MAXORIENTATION ;
	static const double M_T_P_MAXRADIUS  ;
	static const double M_T_P_MAXANGLE  ;
	static const double M_T_P_MAXORIENTATION ;
	static const double M_F_BRUTEFORCE  ;
};
