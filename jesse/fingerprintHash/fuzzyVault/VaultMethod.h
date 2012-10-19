#ifndef __VAULTMETHOD_H__
#define __VAULTMETHOD_H__


//#include "vault.h"

#include "thisLanguageSucks.h"
#include <vector>

class VaultMethod{
	protected:
		vector<float> distro;
		float distroIncrement;
	public:
//		typedef void (VaultMethod::*func)();

		bool INV;

		virtual vector<folded> minutiae2ZZ(const vector<Minutiae> &hTst)=0;
		

		virtual bool compFold(const folded& fold1, const folded& fold2);
		virtual bool comp(const ZZ_p& z1, const ZZ_p& z2);
		virtual float distance(const ZZ_p& z1, const ZZ_p& z2)=0;

//		static virtual bool compFromZero(const folded& fold1, const folded& fold2)=0;
		
		virtual void sorter(const vector<folded> &foldedVector);

		void setDistro(const vector<float>& d);
		vector<float> getDistro();
		void setDistroIncrement(float inc);
		float getDistroIncrement();

  };


class VaultBF : public VaultMethod{
public:
	VaultBF();
	vector<folded> minutiae2ZZ(const vector<Minutiae> &hTst);

	//bool compFold(folded fold1, folded fold2);
	//bool comp(ZZ_p z1, ZZ_p z2);
	float distance(const ZZ_p& z1, const ZZ_p& z2);
//static bool compFromZero(const folded& fold1, const folded& fold2);
};


class VaultTriangles : public VaultMethod{
public:
	VaultTriangles();
	vector<folded> minutiae2ZZ(const vector<Minutiae> &hTst);

	//bool compFold(folded fold1, folded fold2);
	//bool comp(ZZ_p z1, ZZ_p z2);
	float distance(const ZZ_p& z1, const ZZ_p& z2);
};



#endif