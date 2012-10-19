#ifndef __RUN_TEST_H__
#define __RUN_TEST_H__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "Common.h"
#include "Vault.h"
#include "IOUtils.h"

using namespace std;
extern vector<vector<vector<minutia> > > allFingerprints;
extern vector<vector<vector<folded> > > allFingerprintsFolded;
extern vector<vector<vector<vector<folded> > > > allFingerprintsFoldedTransformed;
extern vector<vector<vector<vector<folded> > > > allFingerprintsFoldedRotated;
extern vector<vector<vector<vector<minutia> > > > allFingerprintsRotatedMinu;
extern vector<vector<Vault *> > allVaults;
extern int zz1;
extern int zz2;
//extern vector<float> root;
extern float root[5000000];
vector<vector<Vault * > > dummy();

void makeRoot();
void makeAllVaults(string methodType, float distroIncrement);
void makeAllFolded(string methodType);
void makeAllFoldedTransformed(string methodType);
void makeAllFoldedRotated(string methodType);
vector<float> testGenuines(int fingerprintIndex, int readingNumber, Vault* vault);

vector<float> testImpostors(int fingerprintIndex, Vault* vault);

//@param database is like "FVC2002/DB1"
//@pre 1<=numFingerprints<=100
//optional
float runDatabase(string database, vector<float>* distro, int numFingerprints, float distroIncrement);

double match_tjea(vector<minutia> finger1, vector<minutia> finger2, vector<pair<int,int> > &matchingPairs);

vector<minutia> translate(const vector<minutia> &reference, const vector<minutia> &test, vector<pair<int, int> > &matchingPoints);

#endif
