#ifndef plots_h
#define plots_h
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

#include "dissipativeMap.h"

// Selecting only orbits which stay within y\in [-B,B]
// for M interates.
void plotSystem(DStandardMap &F,int M);
// selecting random trajectories
void plotSystem(DStandardMap &F);

#endif