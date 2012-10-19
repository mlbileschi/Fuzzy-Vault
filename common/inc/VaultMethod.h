#ifndef __VAULTMETHOD_H__
#define __VAULTMETHOD_H__


#include "Common.h"
#include <vector>


class VaultMethod{
	public:

		bool INV;

		virtual std::set<ZZ> minutiae2ZZ(const vector<minutia> &hTst)=0;
		

		virtual bool compFold(const folded& fold1, const folded& fold2);
		virtual bool comp(const ZZ_p& z1, const ZZ_p& z2);
		virtual float distance(const ZZ_p& z1, const ZZ_p& z2)=0;
		virtual void getDefaultDistro(vector<float>* distro, float* distroIncrement)=0;
		virtual string getMethodType()=0;


		virtual std::set<int> minutiae2int(const vector<minutia> &hTst)=0;
		virtual float distance(const int& z1, const int& z2)=0;

  };


class VaultBF : public VaultMethod{
public:
	VaultBF();
	void getDefaultDistro(vector<float>* distro, float* distroIncrement);
	string getMethodType();
	std::set<ZZ> minutiae2ZZ(const vector<minutia> &hTst);

	float distance(const ZZ_p& z1, const ZZ_p& z2);

	std::set<int> minutiae2int(const vector<minutia> &hTst);
	float distance(const int& z1, const int& z2);
};


class VaultTriangles : public VaultMethod{
public:
	VaultTriangles();
	void getDefaultDistro(vector<float>* distro, float* distroIncrement);
	string getMethodType();
	std::set<ZZ> minutiae2ZZ(const vector<minutia> &hTst);

	float distance(const ZZ_p& z1, const ZZ_p& z2);

	std::set<int> minutiae2int(const vector<minutia> &hTst);
	float distance(const int& z1, const int& z2);
};

class VaultLetsTrySomething : public VaultMethod{
public:
	VaultLetsTrySomething();
	void getDefaultDistro(vector<float>* distro, float* distroIncrement);
	string getMethodType();
	std::set<ZZ> minutiae2ZZ(const vector<minutia> &hTst);

	float distance(const ZZ_p& z1, const ZZ_p& z2);

	std::set<int> minutiae2int(const vector<minutia> &hTst);
	float distance(const int& z1, const int& z2);
};


int toFeild(float num, float mean, float sd, int bits);

#endif
