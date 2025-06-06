#ifndef utils_h
#define utils_h
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

inline interval part(interval x, int n, int i) {return x.left()+i*(x.right()-x.left())/n+(x-x.left())/n;}

inline IMatrix toInterval(DMatrix A)
{
	IMatrix B(2,2);
	B[0][0]=A[0][0]; B[0][1]=A[0][1];
	B[1][0]=A[1][0]; B[1][1]=A[1][1];
	return B;
}

inline double toDouble(interval x){return x.mid().leftBound();}

inline double uniformSapmle()
{
	return (double)rand()/(double)RAND_MAX;
}

// This plots the center point of the cartesian product of the 
// intervals, together with the radii of the intervals: rx, ry. 
// The rx, ry can be used in gnuplot to plot boxes.
inline void plot(interval x,interval y,ofstream &file)
{
	double rx=(x.rightBound()-x.leftBound())/2;
	double ry=(y.rightBound()-y.leftBound())/2;
	file << toDouble(x) << " " << toDouble(y) << " " << rx << " " << ry << endl;
}

#endif