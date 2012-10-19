#ifndef IOUTILS_H
#define IOUTILS_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include "../inc/common.h"

using namespace std;

void print(vector<string> v);
void print(vector<float> v);
void print(Minutiae m);
void print(vector<Minutiae> v);

string toString(int number);

vector<string> splitBySpaces(string toSplit);
vector<Minutiae> readFingerprintFromFile(string inputFilename);

void readAllFingerprints(string database);

#endif
