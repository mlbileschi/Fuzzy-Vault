#include "../inc/unitTests.h"
#include "../inc/ProbDistroGA.h"
#include "../inc/evaluatePerformance.h"
#include <iostream>
#include <stdexcept>
#include <vector>

void runUnitTests(){
	int count = 0;
	int pass = 0;

	runVaultUnitTests(&count, &pass);
	runGAUnitTests(&count, &pass);


	cout << "\nResults: " << pass << " of " << count << " test(s) passed\n";
	if(count==pass){cout << "PASSED!!"<<endl;} else { cout << "Your code is bad and you should feel bad"<<endl;}
}
