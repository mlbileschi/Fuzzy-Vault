#ifndef __RUN_TEST_H__
#define __RUN_TEST_H__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "common.h"
#include "vault.h"
#include "IOUtils.h"

using namespace std;
extern vector<vector<vector<Minutiae> > > allFingerprints;
extern vector<vector<vector<folded> > > allFingerprintsFolded;
extern vector<vector<vector<vector<folded> > > > allFingerprintsFoldedTransformed;
extern vector<vector<Vault *> > allVaults;
vector<vector<Vault * > > dummy();
void makeAllVaults(string methodType, float distroIncrement);
void makeAllFolded(string methodType);
void makeAllFoldedTransformed(string methodType);
vector<float> testGenuines(int fingerprintIndex, int readingNumber, Vault* vault);

vector<float> testImpostors(int fingerprintIndex, Vault* vault);

//@param database is like "FVC2002/DB1"
//@pre 1<=numFingerprints<=100
//optional
float runDatabase(string database, vector<float>* distro, int numFingerprints, float distroIncrement);

#endif
