#include "../inc/UnitTests.h"
#include "../inc/ProbDistroGA.h"
#include "../inc/EvaluatePerformance.h"
#include "../inc/RunTest.h"
#include <iostream>
#include <stdexcept>
#include <vector>

//@TODO test for toFeild

bool TestTriangleDistance(int i1, int i2, float expectedDistance){
	cout << "TestTriangleDistance(" << i1 << "," << i2 << "," << expectedDistance << "):\t";
	ZZ_p f1 = to_ZZ_p(i1);
	ZZ_p f2 = to_ZZ_p(i2);
	VaultMethod* method = new VaultTriangles();
	float disto = method->distance(f1, f2);
	if(abs(disto - expectedDistance) <.0001){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed " << disto << endl;
	return false;
}


bool TestBFDistance(int i1, int i2, float expectedDistance){
	cout << "TestBFDistance(" << i1 << "," << i2 << "," << expectedDistance << "):\t\t";
	ZZ_p f1 = to_ZZ_p(i1);
	ZZ_p f2 = to_ZZ_p(i2);
	VaultMethod* method = new VaultBF();
	float disto = method->distance(f1, f2);
	if(abs(disto - expectedDistance) <.0001){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed " << disto << endl;
	return false;
}
		


bool TestToZZTriangles(){
	cout << "TestToZZTriangles:\t\t\t\t";

	vector<minutia> fingerprint;
	minutia next;
	next.x = 215; next.y = 63; next.theta = 263; fingerprint.push_back(next);
	next.x = 198; next.y = 90; next.theta = 181; fingerprint.push_back(next);
	next.x = 185; next.y = 80; next.theta = 43; fingerprint.push_back(next);

	VaultMethod* method = new VaultTriangles();
	std::set<ZZ> zzs = method->minutiae2ZZ(fingerprint);
	if(zzs.count(to_ZZ(47406))){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed ";
	for(std::set<ZZ>::iterator iter=zzs.begin(); iter != zzs.end(); iter++){
		cout << (*iter) << " ";
	}
	cout << endl;	
	return false;
}


bool TestToZZBF(){
	cout << "TestToZZBF:\t\t\t";

	vector<minutia> fingerprint;
	minutia next;
	next.x = 215; next.y = 63; next.theta = 263; fingerprint.push_back(next);
	next.x = 198; next.y = 90; next.theta = 181; fingerprint.push_back(next);
	next.x = 80; next.y = 7; next.theta = 50; fingerprint.push_back(next);
	next.x = 3; next.y = 0; next.theta = 23; fingerprint.push_back(next);

	VaultMethod* method = new VaultBF();
	std::set<ZZ> zzs = method->minutiae2ZZ(fingerprint);
	if(zzs.count(to_ZZ(26746)) && zzs.count(to_ZZ(24759)) && zzs.count(to_ZZ(10242)) && zzs.count(to_ZZ(0))){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed ";
	for(std::set<ZZ>::iterator iter=zzs.begin(); iter != zzs.end(); iter++){
		cout << (*iter) << " ";
	}
	cout << endl;	
	return false;
}



bool Test0(){
	cout << "Test0:\t";
	VaultMethod* method;
//	method = new VaultBF();
	method = new VaultTriangles();
	Vault* vault = new Vault(method);

	//Generate random secret
	vec_ZZ_p secret;
	secret.SetLength(POLYNOMIAL_TERMS);
	for(int i=0; i<POLYNOMIAL_TERMS; i++){
		secret[i] = random_ZZ_p();
	}

	vector<minutia> fingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/FVC2002/DB2/features/27_3.fp");

	vault->lock(fingerprint, secret);
	polyResults result = vault->unlock(fingerprint);
	if(result.score > 10.0){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed"<<endl;
	return false;
}

bool TestSanity0(){
	cout << "TestSanity0:\t\t\t\t";

	VaultMethod* method;
	VaultMethod* method1;
	method = new VaultBF();
	method1 = new VaultTriangles();
	Vault* vault = new Vault(method);
	Vault* vault1 = new Vault(method1);

	//Generate random secret
	vec_ZZ_p secret;
	secret.SetLength(POLYNOMIAL_TERMS);
	for(int i=0; i<POLYNOMIAL_TERMS; i++){
		secret[i] = random_ZZ_p();
	}

	vector<minutia> fingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/FVC2002/DB2/features/27_3.fp");

	vault->lock(fingerprint, secret);
	vault1->lock(fingerprint, secret);
	polyResults result = vault->unlock(fingerprint);
	polyResults result1 = vault1->unlock(fingerprint);
	cout << " -" << result.score << "- ";

	if(result.score > 10.0){
		cout << "passed" << endl;
		return true;
	}
	cout << "failed" << endl;
	return false;
}


bool TestBF(){
	cout << "TestBF:\t\t\t\t";

	VaultMethod* method;
	method = new VaultBF();
	Vault* vault = new Vault(method);

	//Generate random secret
	vec_ZZ_p secret;
	secret.SetLength(POLYNOMIAL_TERMS);
	for(int i=0; i<POLYNOMIAL_TERMS; i++){
		secret[i] = random_ZZ_p();
	}

	vector<minutia> fingerprint;
	vector<minutia> fingerprintTransformed;
	minutia m;
	m.x = 180; m.y = 125; m.theta = 16; fingerprint.push_back(m);
	m.x = 360; m.y = 12; m.theta = 64; fingerprint.push_back(m);
	m.x = 432; m.y = 491; m.theta = 269; fingerprint.push_back(m);
	m.x = 56; m.y = 12; m.theta = 90; fingerprint.push_back(m);
	m.x = 4; m.y = 412; m.theta = 55; fingerprint.push_back(m);

//rotate 6 degrees clockwise about the origin the translate (12,-5)
	m.x = 204; m.y = 101; m.theta = 10; fingerprintTransformed.push_back(m);
	m.x = 371; m.y = -31; m.theta = 58; fingerprintTransformed.push_back(m);
	m.x = 493; m.y = 438; m.theta = 263; fingerprintTransformed.push_back(m);
	m.x = 69; m.y = 1; m.theta = 84; fingerprintTransformed.push_back(m);
	m.x = 59; m.y = 404; m.theta = 49; fingerprintTransformed.push_back(m);


	vault->lock(fingerprint, secret);
	polyResults result = vault->unlock(fingerprintTransformed);
	//cout << " score: " << result.score << "- ";

	if(result.score > 2.9 && result.score < 100.0) {
		cout << "passed" << endl;
		return true;
	}
	cout << "failed: " << result.score <<  endl;
	return false;
}




void runVaultUnitTests(int *count, int *pass){

	makeRoot();
	*pass += TestSanity0(); (*count)++;
	*pass += TestTriangleDistance(10500,10501,0.0); (*count)++;
	*pass += TestTriangleDistance(10500,6098,17.37815); (*count)++;
	*pass += TestTriangleDistance(10500,20226,12.96148); (*count)++;
	*pass += TestTriangleDistance(6098,20227,15.55635); (*count)++;
	*pass += TestToZZTriangles(); (*count)++;
	*pass += TestBFDistance(39908,48188,60.21628); (*count)++;
	*pass += TestBFDistance(39908,3441,52.48809); (*count)++;
	*pass += TestBFDistance(48188,3441,49.56814); (*count)++;
	//*pass += TestToZZBF(); (*count)++;
	*pass += TestBF(); (*count)++;

}
