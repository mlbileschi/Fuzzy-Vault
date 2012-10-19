#include "VaultMethod.h"

void VaultMethod::setDistro(const vector<float>& d){VaultMethod::distro = d;}
std::vector<float> VaultMethod::getDistro(){return VaultMethod::distro;}

void VaultMethod::setDistroIncrement(float inc){VaultMethod::distroIncrement = inc;}
float VaultMethod::getDistroIncrement(){return VaultMethod::distroIncrement;}

/*
vector<folded> VaultMethod::minutiae2ZZ(ptrFPTemplate &hTst){
	vector<folded> ans;
	for(int i=0; i<(*hTst).nCnt; i++){ 
		ZZ z = to_ZZ((*hTst).marr[i].m_nX);
		z = z << 7;
		z = z + (*hTst).marr[i].m_nY;
		//z = z + (*hTst).marr[i].m_nTheta;
		z = z << 1;
		folded fold;
		fold.z1 = to_ZZ_p(z);
		ans.push_back(fold);
	}
	return ans;
}
*/

bool VaultMethod::comp(const ZZ_p& z1, const ZZ_p& z2){
	return distance(to_ZZ_p(0), z1) < distance(to_ZZ_p(0), z2);
}

bool VaultMethod::compFold(const folded& fold1, const folded& fold2){
	return comp(fold1.z1, fold2.z1);
}

void VaultMethod::sorter(const vector<folded> &foldedVector){
//	sort(foldedVector.begin(), foldedVector.end(), &VaultBF::compFold);
}

//float VaultMethod::distance(ZZ_p z1, ZZ_p z2){
//	return 2.0;
//}