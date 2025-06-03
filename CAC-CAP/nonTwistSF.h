#ifndef nonTwistSF_h
#define nonTwistSF_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

// This function validates Theorem 1.3 for a given choice of a
// parameter pair (a,b). The method is essentially the same as 
// the one from conservativeMap.h/cpp. We create a different function
// since due the number of parameter pairs we want to parallelise 
// the computation. In order to do this we pass the maps by reference,
// so that each thread uses differnt map objects.
bool diffusionInNonTwistSF(DMap &f,IMap &F,interval a,interval b,int N_search,int max_n);

#endif