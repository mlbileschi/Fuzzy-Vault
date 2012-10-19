#include "../inc/unitTests.h"
#include "../inc/ProbDistroGA.h"
#include "../inc/evaluatePerformance.h"
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

	vector<Minutiae> fingerprint;
	Minutiae next;
	next.m_nX = 215; next.m_nY = 63; next.m_nTheta = 263; fingerprint.push_back(next);
	next.m_nX = 198; next.m_nY = 90; next.m_nTheta = 181; fingerprint.push_back(next);
	next.m_nX = 185; next.m_nY = 80; next.m_nTheta = 43; fingerprint.push_back(next);

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

	vector<Minutiae> fingerprint;
	Minutiae next;
	next.m_nX = 215; next.m_nY = 63; next.m_nTheta = 263; fingerprint.push_back(next);
	next.m_nX = 198; next.m_nY = 90; next.m_nTheta = 181; fingerprint.push_back(next);
	next.m_nX = 80; next.m_nY = 7; next.m_nTheta = 50; fingerprint.push_back(next);
	next.m_nX = 3; next.m_nY = 0; next.m_nTheta = 23; fingerprint.push_back(next);

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

	vector<Minutiae> fingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/FVC2002/DB2/features/27_3.fp");

	vault->lock(fingerprint, secret);
	polyResults result = vault->unlock(fingerprint);
	if(result.score > 10.0){
		cout << "passed\n";
		return true;
	}
	cout << "failed\n";
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

	vector<Minutiae> fingerprint = readFingerprintFromFile("../../CUBS_FP_DATA/FVC2002/DB2/features/27_3.fp");

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

	vector<Minutiae> fingerprint;
	vector<Minutiae> fingerprintTransformed;
	Minutiae m;
	m.m_nX = 180; m.m_nY = 125; m.m_nTheta = 16; fingerprint.push_back(m);
	m.m_nX = 360; m.m_nY = 12; m.m_nTheta = 64; fingerprint.push_back(m);
	m.m_nX = 432; m.m_nY = 491; m.m_nTheta = 269; fingerprint.push_back(m);
	m.m_nX = 56; m.m_nY = 12; m.m_nTheta = 90; fingerprint.push_back(m);
	m.m_nX = 4; m.m_nY = 412; m.m_nTheta = 55; fingerprint.push_back(m);

//rotate 6 degrees clockwise about the origin the translate (12,-5)
	m.m_nX = 204; m.m_nY = 101; m.m_nTheta = 10; fingerprintTransformed.push_back(m);
	m.m_nX = 371; m.m_nY = -31; m.m_nTheta = 58; fingerprintTransformed.push_back(m);
	m.m_nX = 493; m.m_nY = 438; m.m_nTheta = 263; fingerprintTransformed.push_back(m);
	m.m_nX = 69; m.m_nY = 1; m.m_nTheta = 84; fingerprintTransformed.push_back(m);
	m.m_nX = 59; m.m_nY = 404; m.m_nTheta = 49; fingerprintTransformed.push_back(m);


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
