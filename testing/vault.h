/**
testing purpose: to find the max distance.

*/
#pragma once

#ifndef __VAULT_H__
#define __VAULT_H__


#include "thisLanguageSucks.h"

//Parameters
static float maxDist = -1;
const int NUM_CHAFF_POINTS=200;
const double TRUNK = 1;

const double TOLERANCE = 8; //18
const int KTHRESH = 4;
const int K = 5;


const bool FOLDED = true;
const int POLYNOMIAL_TERMS = 10;

const double THETA_INCREMENT = 5.0;
const double X_INCREMENT = 20.0;
const double Y_INCREMENT = 20.0;

const int THETA_STEPS = 1;
const int X_STEPS = 1;
const int Y_STEPS = 1;

//const double maxDist = 149;

//Fixed
//const int FP_IMAGE_WIDTH=500; 
//const int FP_IMAGE_HEIGHT=500;

const int FP_IMAGE_WIDTH=350; //FVC2002 DB1
const int FP_IMAGE_HEIGHT=350;

const ZZ FEILD_SIZE=to_ZZ(65537);




class Vault{
  protected:
	vector<folded> vault;
	VaultMethod* method;
	ZZ_p function(ZZ_p z, vec_ZZ_p secret);

	public:

		Vault(VaultMethod* m);

		void lock(vector<Minutiae> hRef, vec_ZZ_p secret);

		polyResults unlock(vector<Minutiae> hTst);

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
