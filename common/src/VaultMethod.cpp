#include "../inc/VaultMethod.h"

bool VaultMethod::comp(const ZZ_p& z1, const ZZ_p& z2){
	return distance(to_ZZ_p(0), z1) < distance(to_ZZ_p(0), z2);
}

bool VaultMethod::compFold(const folded& fold1, const folded& fold2){
	return comp(fold1.z1, fold2.z1);
}


//truncates 'num' into a value of 'bits' bits uniformly according a standard distribution given by 'mean' and 'sd'
int toFeild(float num, float mean, float sd, int bits){
	//	mean, mean+.16sd, mean+.32sd, mean+.49sd, mean+.68sd, mean+.89sd, mean+1.15sd, mean+1.53sd, inf .015625
	float scores[] = {-2.15, -1.86, -1.68, -1.53, -1.42, -1.32, -1.23, -1.15, -1.08, -1.01, -.95, -.89, -.83, -.78, -.72,  
			-.68, -.62, -.58, -.53, -.49, -.44, -.40, -.36, -.32, -.28, -.24, -.20, -.16, -.12, -.08, -.04, 
			0, .04, .08, .12, .16, .20, .24, .28, .32, .36, .40, .44, .49, .53, .58, .62, .68, 
			.72, .78, .83, .89, .95, 1.01, 1.08, 1.15, 1.23, 1.32, 1.42, 1.53, 1.68, 1.86, 2.15};
	float zscore = (num - mean)/sd;
	for(int i=0; ((i+1)<<(6-bits))-1<64; i++){
		if(zscore < scores[((i+1)<<(6-bits))-1]){return i;}
	}
	return 63 >> (6-bits);
}	


