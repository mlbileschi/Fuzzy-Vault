//#pragma warning (disable : 4786)

#include "thisLanguageSucks.h"
#include "vault.h"


int Match_Fuzzy_Vault(ptrFPTemplate hRef, ptrFPTemplate hTst, 
							 MatchResult*	pMatchResult, MatchResultEx* pResultEx)
{

	vector<Minutiae> hRefer, hTster;

	for(int i=0; i<hRef->nCnt; i++){
		hRefer.push_back(hRef->marr[i]);
	}
	for(int i=0; i<hTst->nCnt; i++){
		hTster.push_back(hTst->marr[i]);
	}
	
	//cout << (*hRef).nCnt << "  " << " . ";

	ZZ_p::init(FEILD_SIZE);
	pMatchResult->similarity = 0;

	srand(time(NULL));
	SetSeed(to_ZZ(rand()));
	
	//for(int i=0; i<500; i++){
	//	float p = (float)rand()/(float)RAND_MAX;
	//	if(p<0.1){
	//		//cout << p << " ";
	//	}
	////	cout << random_ZZ_p() << " ";
	//}

	VaultMethod* method;
	method = new VaultBF();
//	method = new VaultTriangles();
	Vault* vault = new Vault(method);

	//Generate random secret
	vec_ZZ_p secret;
	secret.SetLength(POLYNOMIAL_TERMS);
	for(int i=0; i<POLYNOMIAL_TERMS; i++){
		secret[i] = random_ZZ_p();
	}



	vault->lock(hRefer, secret);


	// lock vault
	/*
	ofstream file("../../fuckthis.txt", ios::app);
	ofstream filesig("../../fuckthissigma.txt", ios::app);
	for(int i=0; i<hRef->nCnt; i++){
		for(int h=0; h<K; h++){
			file << (*vault)[i].d[h] << "\t";
			filesig << (*vault)[i].sigma[h] << "\t";
		}
		file << "\n";
		filesig << "\n";
	//	file << (*vault)[i].chaff << "\n";
	}

	file.close();
	*/
	//used to test output of the vault
	/* 
	int r = ceil((((double)rand())/((double)RAND_MAX))*100);

	if(r == 86){

	for(int i=0; i<vault->size(); i++){
		cout << "vault x: " << (*vault)[i].x << "\n";
		cout << "vault y: " << (*vault)[i].y << "\n";
		cout << "vault t: " << (*vault)[i].theta << "\n";
		cout << "vault 1: " << (*vault)[i].f1 << "\n";
		cout << "vault 2: " << (*vault)[i].f2 << "\n";
	}
	}
	*/
	/*

	int center_x = FP_IMAGE_WIDTH/2;
	int center_y = FP_IMAGE_HEIGHT/2;


	ptrFPTemplate tester;
	//tester = (ptrFPTemplate) malloc((sizeof(Minutiae)*hTst->nCnt) + 20);
	tester = (ptrFPTemplate) malloc(20000);
	tester->nCnt = hTst->nCnt;

	double max = 0;
	data result;


	//Test various rotations and translations
	int theta_min, theta_max, x_min, x_max, y_min, y_max;

	theta_max = THETA_STEPS/2;
	x_max = X_STEPS/2;
	y_max = Y_STEPS/2;

	if(THETA_STEPS%2 == 1){
		theta_min = -THETA_STEPS/2;
	}else{theta_min = -THETA_STEPS/2 +1;}

	if(X_STEPS%2 == 1){
		x_min = -X_STEPS/2;
	}else{x_min = -X_STEPS/2 +1;}

	if(Y_STEPS%2 == 1){
		y_min = -Y_STEPS/2;
	}else{y_min = -Y_STEPS/2 +1;}

	bool matched_it = false;

	double maxi = -10000000;

	for(int i=theta_min; i<=theta_max; i++){ 
		if(matched_it){break;}
		for(int j=x_min; j<=x_max; j++){
			if(matched_it){break;}
			for(int k=y_min; k<=y_max; k++){
		
			double d_theta = i*THETA_INCREMENT;
			double d_x = j*X_INCREMENT;
			double d_y = k*Y_INCREMENT;

			for(int p=0;p<hTst->nCnt;p++){
				double rot_x = center_x + cos(d_theta/360.0)*(hTst->marr[p].m_nX-center_x) + sin(d_theta/360.0)*(hTst->marr[p].m_nY-center_y);
				double rot_y = center_y - sin(d_theta/360.0)*(hTst->marr[p].m_nX-center_x) + cos(d_theta/360.0)*(hTst->marr[p].m_nY-center_y);
				//cout << "rot: " << (int)rot_x << " " << rot_y << "   ";
				tester->marr[p].m_nX = (int)(rot_x + d_x);
				tester->marr[p].m_nY = (int)(rot_y + d_y);
				tester->marr[p].m_nTheta = hTst->marr[p].m_nTheta + d_theta;
				//cout << "!" << hTst->marr[p].m_nTheta << " " << d_theta << " ";
				//cout << "m: " << tester->marr[p].m_nX << " " << tester->marr[p].m_nY << " " << tester->marr[p].m_nTheta << " ";
				//cout << "m: " << hTst->marr[p].m_nX << " " << hTst->marr[p].m_nY << " " << hTst->marr[p].m_nTheta << " ";
			}

		
			//result = unlockInv(vault, tester); // try to unlock the vault
			result.poly = v;
			result.score = 3;
			if(result.score > maxi){
				maxi=result.score;
			}
	*/
	//

	thisIsTheStructThatIsUsedToReturnMultipleValuesAfterAttemptingToUnlockTheVaultWhichMayExpandInTheFuture result = vault->unlock(hTster);

	if (result.poly == secret){
		cout << "Matched! ";
	}

	//cout << "score: " << result.score << "\n";

	pMatchResult->MatchPerformed = true;
	pMatchResult->similarity = result.score;


	return 0;
}
