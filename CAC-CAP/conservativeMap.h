#ifndef conservativeMap_h
#define conservativeMap_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

void validateChaosInConservativeMaps();
vector<DVector> findShortestTrajectory(DMap &f,double B,int N,int max_n);

#endif