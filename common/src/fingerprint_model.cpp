#include "../inc/fingerprint_model.h"
#include "../inc/random_gaussian.h"

//parameters for numbers of minutia in master fingerprint
const double MEAN_MINUTIA_COUNT=100;
const double STD_MINUTIA_COUNT=20;
//image size 
const int FP_IMAGE_WIDTH=500;
const int FP_IMAGE_HEIGHT=500;
//parameters for the print circle
const double MEAN_PRINT_RADIUS=200;
const double STD_PRINT_RADIUS=20;
const double STD_PRINT_CENTER=20;
//additional affine transform parameters
const double STD_DISPLACEMENT=5; /// was 200 new:25
const double STD_ROTATION_ANGLE=0; /// was 5 new:5
//errors in print minutia positions
const double STD_XY_ERROR=4; // was 7.0
const double STD_THETA_ERROR=3; // was 10
//missing and spurious minutia
const double PROB_MISSING=.2;
const double MEAN_SPURIOUS=0;
const double STD_SPURIOUS=10;

int FingerprintModel::GenerateRandom()
{
	original_finger = vector<minutia>(randn_notrig(MEAN_MINUTIA_COUNT,STD_MINUTIA_COUNT));
//	if(original_finger.size()<0)
//		original_finger.nCnt=0;
//	if(original_finger.nCnt>=CUBS_MAX_MINUTIAE)
//		original_finger.nCnt=CUBS_MAX_MINUTIAE-1;
	int generated=0;
	int n = original_finger.size();
	while(generated < n){
		int x=(double(rand())/double(RAND_MAX))*FP_IMAGE_WIDTH;
		int y=(double(rand())/double(RAND_MAX))*FP_IMAGE_HEIGHT;
		int theta=(double(rand())/double(RAND_MAX))*360;
		bool too_close=false;
		for(int j=0;j<generated;j++){
			double diff_x=original_finger.at(j).x-x;
			double diff_y=original_finger.at(j).y-y;
			if(diff_x*diff_x+diff_y*diff_y<50){
				too_close=true;
				break;
			}
		}
		if(too_close)
			continue;
		original_finger.at(generated).x=x;
		original_finger.at(generated).y=y;
		original_finger.at(generated).theta=theta;
		//original_finger.marr[generated].m_cType='R';
		//original_finger.marr[generated].m_nQuality=1;
		generated++;
	}
	return 0;
}

//Print is a random circle part of the master fingerprint
int FingerprintModel::MakePrint(vector<minutia>& fp)
{
	int center_x, center_y, radius;
	bool done=false;
	while(!done){
		center_x=randn_notrig(FP_IMAGE_WIDTH/2,STD_PRINT_CENTER);
		center_y=randn_notrig(FP_IMAGE_HEIGHT/2,STD_PRINT_CENTER);
		//cout << center_x << "|";
		radius=randn_notrig(MEAN_PRINT_RADIUS,STD_PRINT_RADIUS);
		//cout << "|" << radius << "|";
		if(center_x-radius>0 && center_x+radius<FP_IMAGE_WIDTH
			&& center_y-radius>0 && center_y+radius<FP_IMAGE_HEIGHT)
			done=true;
	}
	int copied=0;
	for(int i=0;i<original_finger.size();i++){
		//cout << "3";
		int diff_x=original_finger.at(i).x-center_x;
		int diff_y=original_finger.at(i).y-center_y;
		//cout << "    " << original_finger.at(i).y << " " << center_x << " " << diff_x*diff_x <<" "<< diff_y*diff_y << " " << radius*radius << endl;
		if(diff_x*diff_x + diff_y*diff_y<radius*radius){
			//cout << "2";
			if((double(rand())/double(RAND_MAX))<PROB_MISSING) 
				continue; //this minutia, though inside the circle, will be missing
			fp.push_back(original_finger.at(i));
			//cout << "1";
//			fp.at(copied).m_nFeatureNo=copied;
			copied++;
		}
	}

	//cout << " " <<  fp.size() << "  " ;
	//cout << copied << "  " ;
	//add spurious minutia
	int num_spurious;
	while(num_spurious=randn_notrig(MEAN_SPURIOUS,STD_SPURIOUS)<0);
	if(copied+num_spurious>=CUBS_MAX_MINUTIAE)
		num_spurious=CUBS_MAX_MINUTIAE-1-copied;
	while(num_spurious>0){
		int x=(double(rand())/double(RAND_MAX))*FP_IMAGE_WIDTH;
		int y=(double(rand())/double(RAND_MAX))*FP_IMAGE_HEIGHT;
		int theta=(double(rand())/double(RAND_MAX))*360;
		double diff_x=x-center_x;
		double diff_y=y-center_y;
		if(diff_x*diff_x+diff_y*diff_y>=radius*radius) //not inside print circle
			continue;
		bool too_close=false;
		for(int j=0;j<copied;j++){
			diff_x=fp.at(j).x-x;
			diff_y=fp.at(j).y-y;
			if(diff_x*diff_x+diff_y*diff_y<50){
				too_close=true;
				break;
			}
		}
		if(too_close)
			continue;
		fp.at(copied).x=x;
		fp.at(copied).y=y;
		fp.at(copied).theta=theta;
//		fp.marr[copied].m_cType='R';
//		fp.marr[copied].m_nQuality=1;
		copied++;
		num_spurious--;
	}

	//introduce errors in minutia positions and direction
	for(int i=0;i<fp.size();i++){
		fp.at(i).x+=randn_notrig(0,STD_XY_ERROR);
		fp.at(i).y+=randn_notrig(0,STD_XY_ERROR);
		fp.at(i).theta+=randn_notrig(0,STD_THETA_ERROR);
	}

	
	//apply global affine transformation
	double d_theta=randn_notrig(0,STD_ROTATION_ANGLE);
	double d_x=randn_notrig(0,STD_DISPLACEMENT);
	double d_y=randn_notrig(0,STD_DISPLACEMENT);
	for(int i=0;i<fp.size();i++){
		double rot_x=center_x+cos(d_theta/360)*(fp.at(i).x-center_x)+sin(d_theta/360)*(fp.at(i).y-center_y);
		double rot_y=center_y-sin(d_theta/360)*(fp.at(i).x-center_x)+cos(d_theta/360)*(fp.at(i).y-center_y);
		fp.at(i).x=rot_x+d_x;
		fp.at(i).y=rot_y+d_y;
		fp.at(i).theta=fp.at(i).theta+d_theta;
	}
	return 0;
}
