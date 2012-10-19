#ifndef EVALUATE_PERFORMANCE_H
#define EVALUATE_PERFORMANCE_H
#include <vector>
#include "Common.h"

typedef struct
{
	float EER;
	float TER;
	float FAR01;
	float FAR001;
	float FAR0001;
	float FAR0;
} ROC;

int is_this_working();

float evaluatePerformance(vector<float> impostors, vector<float> genuines);

#endif // EVALUATE_PERFORMANCE_H



