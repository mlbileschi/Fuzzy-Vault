#include <math.h>
#include "../inc/matching.h"
#include "../inc/SFeature.h"
/*
  File : polar_util.cpp
*/

void cart2polar(
  int cx,       // The current point x coordinate.
  int cy,       // The current point y coordinate.
  int rx,       // The reference point x coordinate.
  int ry,       // The reference point y coordinate.
  double* r,    // [Output] radius
  double* theta // [Output] angle (0 <= angle < 360)
  ){
  double x, y;
  x = (double)(cx - rx);
  y = (double)(cy - ry);
	*r = 0.0;
	*theta = 0.0;
	*r = sqrt(x * x + y * y);
	if(x > 0.0){
		if(y > 0.0){
			*theta = atan(y/x) * 180 / PI;
		}else if(y < 0.0){
      *theta = 360 - (atan(ABSI(y/x)) * 180 / PI);
		}else{
			*theta = 0.0;
		}
	}else if(x < 0.0){
		if(y > 0.0){
      *theta = 180.0 - (atan(ABSI(y/x)) * 180 / PI);
		}else if(y < 0.0){
			*theta = atan(ABSI(y/x)) * 180 / PI + 180;
		}else{
			*theta = 180.0;
		}
	}else{
	  if(y > 0.0)
		  *theta = 90.0;
		else if (y < 0.0)
		  *theta = 270.0;
		else
		  *theta = 0.0;
	}
}

/* return the absolute angle difference between 0 and 180 */
double abs_angle_diff(double x, double y){
  double result;
  result = ABSD(x - y);
  if(result > 180.0) return 360 - result;
  else return result;
}
