#ifndef IOUTILS_H
#define IOUTILS_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include "../inc/Common.h"

using namespace std;

string concatenateNotFirst(int len, char* arr[]);

void print(vector<string> v);
void print(vector<float> v);
void print(minutia m);
void print(vector<minutia> v);

string toString(int number);

vector<string> splitBySpaces(string toSplit);
vector<minutia> readFingerprintFromFile(string inputFilename);

void readAllFingerprints(string database);
void center(vector<minutia>& print);
void makeAllRandomFingerprints();

#endif
