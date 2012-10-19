/** 
@author: Max Bileschi
email: mlbileschi@gmail.com

This file treats a vector<float> as a template for a genome for a population, and
runs a genetic algorithm on that population. The standard utilities are implemented,
including an evaluation of fitness, calculations to produce a new generation, 
*/
#include "../inc/ProbDistroGA.h"
#include "../inc/RunTest.h"

#include <fstream>
#include <iterator>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <float.h>
using namespace std;



/*	runs a genetic algorithm on a population of size of @param population size
	with a genome (which is a vector) of size @param genomeSize. The algorithm
	will run for @param iterations steps. Each iteration, it will report the 
	percent completed, the best for that generation, and the variance of the 
	population, which measures "convergence" in some rough sense. 
	@param database is the fingerprint database we're looking at.
*/
void beginGenetic(int populationSize, int genomeSize, float distroIncrement, int iterations, 
				string database, int numFingerprints, string call)
{
	srand(time(NULL));
	vector<vector< float> > generation = newGeneration(populationSize, genomeSize);
	cout<<"BEGIN GENETIC ALGORITHM"<<endl;
	vector<float> best(genomeSize);
	float fitnessOfBest = 0;

	for(int i = 0; i<iterations; i++)
	{
		//spawn a new generation: reproduce, mutate, crossover
		//GenProps is in file ProbDistroGA.h
		GenProps genProps = nextGeneration(generation, populationSize, genomeSize, 
							numFingerprints, distroIncrement, database);
		generation = genProps.generation;

		//remember best fitness so far, and its genome
		if(genProps.fitnessOfBest>fitnessOfBest)
		{
			fitnessOfBest = genProps.fitnessOfBest;
			best = genProps.best;
		}

		//updates, instead of printing a new line (due to \r)
		//prints fitness (we actually want to look at 1/fitness because EER 
		// is better when it's low)
		cout<<"\r"<<"fitness of best = "<<1.0/fitnessOfBest<<"\t"<<
					"convergence : "<<genProps.convergence<<
					"\t  "<<(int)(((float)i/iterations)*100)<<"\% finished";
		cout.flush();
		writeGenomeToFile("../results/"+call+".txt", best, 1/fitnessOfBest);
	}

	cout<<endl<<"END GENETIC ALGORITHM"<<endl;
	cout<<endl<<"best was: "<<endl;
	printDistro(best);
	cout<<"with score : "<<1.0/fitnessOfBest<<endl; //again, low EER is better
	writeGenomeToFile("../results/"+call+"_FINISHED.txt", best, 1/fitnessOfBest);
}



//just pretty-prints @param d
void printDistro(vector<float> d)
{
	cout<<"(";
	for(int i = 0; i<d.size(); i++)
	{
		cout<<d.at(i)<<",";
	}
	cout<<")"<<endl;
}

//just pretty-prints @param g
void printGeneration(vector< vector <float> > g)
{
	for(unsigned int i = 0; i<g.size(); i++)
	{
		printDistro(g.at(i));
	}
	cout<<endl<<endl;
}

/*
Creates a new vector of size @param size of floats 0<=num<=1
that are decreasing
*/
vector<float> newDistro(int size)
{
	vector<float> toret;

	for(int i=0; i<6; i++){ //    13     <-----  |||*** Adjust distro here ***|||
		toret.push_back(1.0);
	}
	toret.push_back(0.0);

/*
	for(unsigned i = 0; i<size; i++)
	{
		toret.push_back((float)rand()/((float)RAND_MAX));
	}

	sort (toret.begin(), toret.end(),std::greater<float>()); //sort descending
*/	
	return toret;
}

bool writeGenomeToFile(string file, vector<float> genome, float fitness)
{
	ofstream vectorFile(file.c_str(), ios::out);
	vector<float>::iterator i;
	for (i = genome.begin(); i < genome.end(); ++i) {
		vectorFile << *i;
		vectorFile << ", ";
	}
	vectorFile<<endl<<"fitness: "<<fitness<<endl;

	return true;
}

/*
Creates a new vector of size @param size of floats 0<=num<=1
that are decreasing
*/
vector<float> newThresholdDistro(int size)
{
	vector<float> toret;
	int numOnes = (int)(ceil(size*(float)rand()/((float)RAND_MAX)))+1;

	for(unsigned i = 0; i<numOnes; i++)
	{
		toret.push_back(1.0);
	}
	for(unsigned i = 0; i<size-numOnes; i++)
	{
		toret.push_back(0.0);
	}

	return toret;
}

/*
Creates @param number new vectors of size @param size, by 
calling newDistro @param number times.
*/
vector< vector <float> > newGeneration(int number, int size)
{
	vector< vector<float> > toret;
	for(unsigned int i = 0; i<number; i++)
	{
		toret.push_back(newDistro(size));
	}
	return toret;
}

/*
just a utility function to take the average over all entries of
@param sample
*/
float mean(vector<float> sample)
{
    float sum = 0;
    for(int i = 0; i < sample.size(); i++)
        sum += sample.at(i);
    return (sum / sample.size());
}

/*
just a utility function to compute the variance of a sample
*/
float variance(vector<float> sample)
{
    float avg = mean(sample);

    float temp = 0;
    for(int i = 0; i < sample.size(); i++)
    {
         temp += (sample.at(i) - avg) * (sample.at(i) - avg) ;
    }
    return temp / sample.size();
}

/*
measures the variance of each coordinate of the genome of each in a generation,
and then takes the mean. It's not particularly accurate, but unfortunately, 
there's not even a consensus in the literature on how to measure
convergence
*/
float convergenceOfGeneration(vector<vector<float> > generation)
{
	vector<float> variances;

	//assumes all genomes are same size
	//for each entry in a genome, visit every genome in generation
	for(int i = 0; i<generation.at(0).size(); i++)
	{
		vector<float> coordinate;  //a coord of each animal in a generation, 
									//e.g. generation[0][1], generation[1][1], ...
		for(int j = 0; j<generation.size(); j++)
		{
			coordinate.push_back(generation.at(j).at(i));
		}
		variances.push_back(variance(coordinate));
	}

	return 1/mean(variances);
}

/*
Uses 1/EER (equal error rate) of a fingerprint reading to determine
fitness. runDatabase is in file RunTest.cpp
*/
float evalFitness(vector<float> genome, int numFingerprints, float distroIncrement, string database)
{
	return 1.0/runDatabase(database, &genome, numFingerprints, distroIncrement);
	//return 1;
}

/*
Takes two genomes, and "mixes up" 1/10th of them. (can mix up same spot twice)
Then it returns a sorted version.
*/
vector<float> crossover(vector<float> genome1,vector<float> genome2)
{
	vector<int> toMutate;
	for(unsigned int i = 0; i<genome1.size()/10; i++)
	{
		toMutate.push_back(rand()%genome1.size()); //rand number between 0 and 
	}
	for(int i = 0; i<toMutate.size(); i++)
	{
		genome1.at(toMutate.at(i)) = (genome1.at(toMutate.at(i))
											+ genome2.at(toMutate.at(i)))/2;
	}
	sort (genome1.begin(), genome1.end(),std::greater<float>()); //sort descending
	return genome1;
}

/*
takes in a genome, and randomizes 1/10th of it (can affect same spot twice)
returns sorted
*/
vector<float> mutate(vector<float> genome)
{
	vector<float> mutated=genome;
	for(int i = 0; i<genome.size()/10; i++)
	{
		mutated.at(rand()%genome.size()) = ((float)rand()/(float)(RAND_MAX)); //rand number between 0 and 99
	}
	sort (mutated.begin(), mutated.end(),std::greater<float>()); //sort descending
	return mutated;
}

/*returns a new generation, best one from the previous generation, the fitness of the best, 
and the relative convergence of solution
NOTE: low EER is a good thing, and we are using 1/EER as a measure of fitness
so current metric is 1/eer for determining how many progeny
*/
GenProps nextGeneration(vector< vector <float> > g, int populationSize, int genomeSize, 
								int numFingerprints, float distroIncrement, string database)
{
	GenProps genProps;
	float totalSuccess = 0;
	vector<float> fitnesses;


	for(int i = 0; i<g.size(); i++)
	{
		vector<float> current = g.at(i);
		fitnesses.push_back(evalFitness(current, numFingerprints, distroIncrement, database));
		totalSuccess += fitnesses.back();
	}

	float maxFitness = fitnesses.at(0);
	int idxOfMaxFitness = 0;
	//find best fitness
	for(int i = 0; i<fitnesses.size()-1; i++)
	{
		if(fitnesses.at(i)>maxFitness)
		{
			maxFitness = fitnesses.at(i);
			idxOfMaxFitness = i;
		}
	}

	vector<float> best = g.at(idxOfMaxFitness);

	vector<vector<float> > newGen;
	for(int i = 0; i<g.size()-1; i++)
	{
		vector<float> current = g.at(i);

		//why + (constant) ? it affects how quickly the solution converges
		//i.e. fitnesses.at(i)/totalSuccess*g.size() represents the total 
		//percentage of the fitness of an entire population had by 
		//fitnesses.at(i). 
		for(int k = 0; k<floor(fitnesses.at(i)/totalSuccess*g.size())+0.0; k++)
		{
			if(rand()%5==0)//one in a twenty times it mutates
				newGen.push_back(mutate(current)); 
			else if (rand()%5==1)//one in twenty times it crosses over
				newGen.push_back(crossover(g.at(rand()%g.size()), current));
			else
				newGen.push_back(current);
		}
	}

	//because we add a little to floor k in loop above
	while(newGen.size()>populationSize)
	{
		newGen.pop_back();
	}
	
	while(newGen.size()<populationSize)
	{
		newGen.push_back(newDistro(genomeSize));
	}

	genProps.convergence = convergenceOfGeneration(newGen);
	genProps.generation = newGen;
	genProps.best = best;
	genProps.fitnessOfBest = maxFitness;
	
	return genProps;
}
