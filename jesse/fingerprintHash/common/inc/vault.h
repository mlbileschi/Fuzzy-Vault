#pragma once

#ifndef __VAULT_H__
#define __VAULT_H__


#include "common.h"

//Parameters

const int NUM_CHAFF_POINTS=200;
const double TRUNK = 1;

const double TOLERANCE = 8; //18
const int KTHRESH = 4;
const int K = 5;


const bool FOLDED = true;
const int POLYNOMIAL_TERMS = 10;

const double THETA_INCREMENT = 15.0;
const double X_INCREMENT = 8.0;
const double Y_INCREMENT = 8.0;

const int THETA_STEPS = 3;
const int X_STEPS = 5;
const int Y_STEPS = 5;

//const double maxDist = 149;

//Fixed
//const int FP_IMAGE_WIDTH=500; 
//const int FP_IMAGE_HEIGHT=500;

const int FP_IMAGE_WIDTH=350; //FVC2002 DB1
const int FP_IMAGE_HEIGHT=350;

const ZZ FEILD_SIZE=to_ZZ(65537);




class Vault{
	private:
		float maxDistro;

		float distroIncrement;
		vector<float> distro;

	protected:

		vector<folded> vault;
		VaultMethod* method;
		ZZ_p function(ZZ_p z, vec_ZZ_p secret);

	public:


		Vault(VaultMethod* m);

		void lock(vector<Minutiae> hRef, vec_ZZ_p secret);

		polyResults unlock(const vector<vector<folded> >& hTst);
		polyResults unlock(const vector<Minutiae>& hTst);
		polyResults unlock(const vector<folded>& hTst);
		void keyInLock(polyResults* data, const std::set<int>& indexPoints);

		std::set<int> compareExtract(polyResults* data, const vector<Minutiae>& hTst);
		std::set<int> compareExtract(polyResults* data, const vector<folded>& hTst);
		std::set<int> compareExtractBF(polyResults* data, const vector<vector<folded> >& translations);
		std::set<int> compareExtractBF(polyResults* data, const vector<folded>& hTst);
		std::set<int> compareExtractBF(polyResults* data, const vector<Minutiae>& hTst);

		void setDistro(const vector<float>& d);
		vector<float> getDistro();
		void setDistroIncrement(const float& inc);
		float getDistroIncrement();
		
		void computeMaxDistro();

		bool operator()(const folded& f1, const folded& f2);

}; // end class

//this class was added to get c++ to do the simple task of sorting properly
class Compare_functor {
	public:
		Compare_functor(VaultMethod* g_) : g(g_) {}
		bool operator()(const folded& a, const folded& b) {
			return g->compFold(a, b);
		}
	private:
		VaultMethod* g;
};
	 



#endif //__VAULT_H__
