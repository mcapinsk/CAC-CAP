#ifndef conservativeMap_h
#define conservativeMap_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

// This routine provides a computer assisted proof of Theorem 1.2 
// from the paper. 
void validateChaosInConservativeMaps();

// This routine finds a trajectory which goes from level y=-B to
// above y=B. This routine is not rigorous. Such trajectory plays the
// role of an initial guess, which is then validated by means of the 
// parallel shooting Krawczyk method from validatedTrajectory.h/cpp
vector<DVector> findShortestTrajectoryUp(DMap &f,double B,int N,int max_n);

#endif