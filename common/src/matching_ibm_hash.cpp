// Author:  Maulik Dave (CSE 666 ) Spring 2012.
// Based on the paper "Anonymous and Revocable Fingerprint Recognition" IBM TJ Watson Research Center (Farooq 2007)
// This algorithm introduces binary strings representation of fingerprints.

#pragma warning (disable : 4786)

#include "common.h"
#include "matching.h"
#include "global_parameters.h"
#include "SFeature.h"
#include "logfile.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>

using namespace std;

// Map data structure will be used for histogram - index and counter.
map<int,int> Ref_hist;   // ??? How to initialize counter value to be zero. 
map<int,int> Test_hist;

//--------------------------------------------------------------------------------------------
//cart_dist between two points (x1, y1) and (x2, y2)
//--------------------------------------------------------------------------------------------
double cart_dist(int x1, int y1, int x2, int y2){
    return sqrt((double)((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

//---------------------------------------------------------------------------------------------
// conv_to_binary_string takes 7 invariants (feature vector from a triangle) and Converts
// it to 24 bits string as explained in code.
// Calcuate this 24 bit number for al the triples and histogram of length n = 2^24 for all the triangles.
//----------------------------------------------------------------------------------------------
int conv_to_binary_string(double s1, double s2, double s3, double a1, double a2, double a3, double h) {
	
	int index  = 0;

	// Quantization - side (S1,S2,S3) 4 bits, angle (A1,A2,A3) 3 bits and triangle height (H) 3 bits.
	// A1 - A3 range 0-359 Degree. Quantizing them by 3 bits, 000-111 (0-7).
	// 360/8 = 45
	int A1 = static_cast<int>(a1/45);
    int A2 = static_cast<int>(a2/45);
    int A3 = static_cast<int>(a3/45);
	
	//S1-S3 Range 0-500, Quantizing them by 4 bits, 0000-1111 (0-15)
	// 500/16 = 31.25
	int S1 = static_cast<int>(s1/31.25);
	int S2 = static_cast<int>(s2/31.25);
	int S3 = static_cast<int>(s3/31.25);

	//H Range 0-500, Quantizing them by 3 buts, 000-111 (0-7)
	// 500/8 = 62.5
	int H  = static_cast<int>(h/62.5);

    // Append the quantized invariants A1A2A3S1S2S3H is 24 bit string.
    /*index is 32 bit integer....thinking like it is 32 bit array:
    index[32:24] = 0;
    index[23:21] = binary value of A1
    index[20:18] = binary value of A2
    index[17:15] = binary value of A3
    index[14:11] = binary value of S1
    index[10:7] = binary value of S2
    index[6:3] = binary value of S3
    index[2:0] = binary value of H*/

	index = ((A1 << 21) | (A2 << 18) | (A3 << 15) | (S1 << 11) | (S2 << 7) | (S1 << 3) | (H << 0));
	return index;

    // My first attempt to quantize with Masking.
	//index = static_cast<int>(a1) & 0x00E00000;
	//index = index || (static_cast<int>(a2) & 0x00000000);
}

//-------------------------------------------------------------------
// Check if given three points it makes valid triangle.
//-------------------------------------------------------------------
int make_triangle(int x1, int y1, int t1, int x2, int y2, int t2, int x3, int y3, int t3) {
 		
		double s1 = cart_dist(x1,y1,x2,y2);
		double s2 = cart_dist(x1,y1,x3,y3);
		double s3 = cart_dist(x2,y2,x3,y3);
        double s = (s1+s2+s3)/2;
		double A = sqrt(s*(s-s1)*(s-s2)*(s-s3)); // Heron's Formula to find the Area of a triagle.
        
		double max_side;
		if (s1 > s2)
			max_side = s1;
		else
			max_side = s2;
		if (s3 > max_side)
			max_side = s3;

		double h = 2*A/max_side;  

		double rad = (180.0 /3.14159265);

		double alpha =  acos( (s2*s2 + s1*s1 - s3*s3)/(2*s1*s2)) * rad;
        double beta  =  acos( (s1*s1 + s3*s3 - s2*s2)/(2*s1*s3)) * rad;
        double gamma =  acos( (s2*s2 + s3*s3 - s1*s1)/(2*s2*s3)) * rad;

		double a1 = t1 - atan2(double (y2-y1), double (x2-x1)) * rad ;
		double a2 = t2 - atan2(double (y3-y2), double (x3-x2)) * rad ;
		double a3 = t3 - atan2(double (y1-y3), double (x1-x3)) * rad ;

		// making sure angles are in the range of 0 - 359 degrees.
		if (a1 < 0)
			a1 = a1 + 360;
		if (a1 >= 360)
			a1 = a1 - 360;

		if (a2 < 0)
			a2 = a2 + 360;
		if (a2 >= 360)
			a2 = a2 - 360;

		if (a3 < 0)
			a3 = a3 + 360;
		if (a3 >= 360)
			a3 = a3 - 360;

		// Commenting out as it may be redundant to check this.
		//if ((x1 == x2) && (x1 == x3) ) {
		//	cout << "In valid Triangle - three points in same line." << "\n";
		//    return -1;
		//}
		//if ((y1 == y2) && (y1 == y3) ) {
		//	cout << "In valid Triangle - three points in same line." << "\n";
		//    return -1;
		//}

		// Tringle internal angles add up to 180, doing this way to avoid double accuracy issue.
		double angle_sum = alpha + beta + gamma - 180.0;
        if (angle_sum  < 0)
			angle_sum = -angle_sum;
		if ( angle_sum > 0.001)  {
			//cout << "In valid Triangle - Angles don't add up to 180." << "\n";
		    return -1;
		}
		// Triangle sum of the two sides can't be greater than third.
		if ( (s1 > s2+s3) || (s2 > s1+s3) || (s3 > s1+s2) ) {
			//cout << "In valid Triangle - sum of two sides greater than third side." << "\n";
		    return -1;
		}
		int str = conv_to_binary_string(s1,s2,s3,a1,a2,a3, h);
		return str;

		/*cout << "Triangle: " << "\n";
		cout << "Pt1(x,y,theta): (" << x1 << " " << y1 << " " << t1 << ")\n";
		cout << "Pt2(x,y,theta): (" << x2 << " " << y2 << " " << t2 << ")\n";
		cout << "Pt3(x,y,theta): (" << x3 << " " << y3 << " " << t3 << ")\n";
		cout << "alpha(internal angle 1): " << alpha << "\n";
		cout << "beta(internal angle 2): " << beta << "\n";
		cout << "gamma(internal angle 3): " << gamma << "\n\n";
		cout << "s1: " << s1 << "\n";
		cout << "s2: " << s2 << "\n";
		cout << "s3: " << s3 << "\n";
		cout << "a1: " << a1 << "\n";
		cout << "a2: " << a2 << "\n";
		cout << "a3: " << a3 << "\n";
		cout << "h: " << h << "\n";
		cout << "string: " << str <<"\n\n";*/
}
//-------------------------------------------------------------------
// Find Similarity between Reference and Test Fingerprint
// May be Binarize here to avoid duplication of travesring ?????
//-------------------------------------------------------------------
double find_similarity() {

	map<int,int>::iterator it1;
	map<int,int>::iterator it2;
	int Ref_one = 0;
	int Test_one = 0;
	int common_one = 0;
	double sim = 0.0;

	it2 = Test_hist.begin();
	for (it1 = Ref_hist.begin(); it1 != Ref_hist.end(); it1++) {	

		if ((*it1).second == 1) { 

			Ref_one++;
			// Do Ref and Test have "1" at same key positions?
			it2 = Test_hist.find((*it1).first); 
			if ((it2 != Test_hist.end()) && ((*it2).second == 1))
				common_one++;
			
		} // outer if
	} // for iteration
   
	for (it2 = Test_hist.begin(); it2 != Test_hist.end(); it2++) {	
		if ((*it2).second == 1) 
			Test_one++;
	}

	cout << "\n" << "Number of 1 at same position: " << common_one << "\n";
	cout << "Total number of 1 in Reference: " << Ref_one << "\n";
	cout << "Total number of 1 in Test: " << Test_one << "\n";

	if (Ref_one == 0)
		return 0.0;
	if (Test_one == 0)
		return 0.0;
	
	double norm = (double) Ref_one*Test_one;
	sim = common_one / sqrt(norm);
	return sim;
	
} // find_similarity

//-------------------------------------------------------------------
// Finger print matching between hRef (subject) and hTst (test)
//-------------------------------------------------------------------
int CUBS_Match_IBM_Hash(ptrFPTemplate hRef, ptrFPTemplate hTst, 
							   MatchResult*	pMatchResult, MatchResultEx* pResultEx)
{
    Minutiae* pMRef;    // pointer to Minutia structure for Reference fingerprint input
    int       nCntRef;  // number of valid Minutia points for Reference fingerprint input

	Minutiae* pMTest;   // pointer to Minutia structure for Test fingerprint input
    int       nCntTest; // number of valid Minutia points for Reference fingerprint input
	int i, j, k;        // for loop indexes
	int n;              // number of minutias
	map<int,int>::iterator it; // iterator to add triangle index and histogram.

	// Basic error checking adequate minutia for Reference and Test input.
	if(hRef==NULL || hRef->nCnt < 2 || hRef->nCnt > CUBS_MAX_MINUTIAE )
		return ERR_MISSING_MINUTIAE;
	if(hTst==NULL || hTst->nCnt < 2 || hTst->nCnt > CUBS_MAX_MINUTIAE )
		return ERR_MISSING_MINUTIAE;

	// Access array of Minutia from Reference fingerprint input
	pMRef   = (Minutiae*)hRef->marr;
    nCntRef = hRef->nCnt;
    // Access array of Minutia from Test fingerprint input
    pMTest   = (Minutiae*)hTst->marr;
    nCntTest = hTst->nCnt;

    // Creating valid triangles out of these minutia points for Reference
	n = nCntRef;   
	for(i=0; i< n; i++) {
		for (j = i+1; j < n; j++) {
			for (k = j+1; k < n; k++) {
		
				int x1 = pMRef[i].m_nX;
				int y1 = pMRef[i].m_nY;
				int t1 = pMRef[i].m_nTheta;
			    int x2 = pMRef[j].m_nX;
				int y2 = pMRef[j].m_nY;
				int t2 = pMRef[j].m_nTheta;	
				int x3 = pMRef[k].m_nX;
				int y3 = pMRef[k].m_nY;
				int t3 = pMRef[k].m_nTheta;

				//cout << "Reference Triangle" << "\n";
		        int str = make_triangle(x1,y1,t1,x2,y2,t2,x3,y3,t3);
				// if str <=0 not a valid triangle.
				if (str > 0){ 
					it = Ref_hist.find(str);
					if(it == Ref_hist.end())
						Ref_hist.insert(pair<int,int>(str,0));
					else
						 Ref_hist[str] = Ref_hist[str] + 1; 
				} // if/else
			}
		}
	} // End of for loop for creating triangles for Reference 

    // Creating valid triangles out of these minutia points for Test.
	n = nCntTest;                 
	for(i=0; i< n; i++) {
		for (j = i+1; j < n; j++) {
			for (k = j+1; k < n; k++) {
		
				int x1 = pMTest[i].m_nX;
				int y1 = pMTest[i].m_nY;
				int t1 = pMTest[i].m_nTheta;
			    int x2 = pMTest[j].m_nX;
				int y2 = pMTest[j].m_nY;
				int t2 = pMTest[j].m_nTheta;	
				int x3 = pMTest[k].m_nX;
				int y3 = pMTest[k].m_nY;
				int t3 = pMTest[k].m_nTheta;

				//cout << "Test Triangle" << "\n";
		        int str = make_triangle(x1,y1,t1,x2,y2,t2,x3,y3,t3);
				// if str <= 0, not a valid triangle.
				if (str > 0){ 
					it = Test_hist.find(str);
					if(it == Test_hist.end())
						Test_hist.insert(pair<int,int>(str,0));
					else
						 Test_hist[str] = Ref_hist[str] + 1; 
				} // if/else 
			}
		}
	} // End of for loop for creating triangles for Test

	// Check Reference/Test Histogram values.	
	/*cout << " Reference Histogram " << "\n";
	for (it = Ref_hist.begin(); it != Ref_hist.end(); it++) {	
		cout << (*it).first << " => " << (*it).second << endl;
	}
    cout << " Test Histogram " << "\n";
	for (it = Test_hist.begin(); it != Test_hist.end(); it++) {	
		cout << (*it).first << " => " << (*it).second << endl;
	}*/

    // Binarize this histogram as for all i, if Fi = 0 if Fi not 1, else Fi = 1
	//cout << " Reference Histogram Binarized" << "\n";
	for (it = Ref_hist.begin(); it != Ref_hist.end(); it++) {	
		if ((*it).second != 1 )
			(*it).second = 0;
		else 
			(*it).second = 1;
		//cout << (*it).first << " => " << (*it).second << endl;
	}
    //cout << " Test Histogram Binarized " << "\n";
	for (it = Test_hist.begin(); it != Test_hist.end(); it++) {	
		if ((*it).second != 1 )
			(*it).second = 0;
		else
			(*it).second = 1;
		//cout << (*it).first << " => " << (*it).second << endl;
	}

	// Output, Matching similarity.
	pMatchResult->MatchPerformed = true;
	pMatchResult->similarity = 1 - find_similarity();
	cout << "Score: " << pMatchResult->similarity << "\n";

    return ERR_NO_ERROR;
}