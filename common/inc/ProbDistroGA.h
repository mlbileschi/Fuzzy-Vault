#ifndef __PROB_DISTRO_GA_H__
#define __PROB_DISTRO_GA_H__


#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
using namespace std;

class GenProps
{
	public: vector<vector<float> > generation;
		vector<float> best;
		float fitnessOfBest;
		float convergence;
};

float mean(vector<float> sample);

float variance(vector<float> sample);

vector<float> newDistro(int genomeSize);

float convergenceOfGeneration(vector<vector<float> > generation);

// @param database looks like "FVC2002/DB1"
// call is what the whole program was run with, and will help to discern separate results
// from different runs of the program. a file will be written to named "call"
void beginGenetic(int populationSize, int genomeSize, float distroIncrement, int iterations, 
				string database, int numFingerprints, string call);

void printGeneration(vector<vector<float> >);

void printDistro(vector<float> );

vector<float> newRandomDistro(int size);

vector<float> newThresholdDistro(int size);

vector< vector <float> > newGeneration(int number, int size);

bool writeGenomeToFile(string file, vector<float> genome, float fitness);

void printDistro(vector<float> d);

float evalFitness(vector<float> probDistro, int numFingerprints, float distroIncrement, string database);

vector<float> crossover(vector<float> probDistro1,vector<float> probDistro2);

vector<float> mutate(vector<float> probDistro);

GenProps nextGeneration(vector< vector <float> >, int populationSize, 
									int genomeSize, int numFingerprints, float distroIncrement, string database);

#endif //__PROB_DISTRO_GA_H__
