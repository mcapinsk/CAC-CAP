#ifndef nonTwistSF_h
#define nonTwistSF_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

interval diffusionAreaNonTwistSF(interval A,interval B,int N_search,int max_n);
bool diffusionInNonTwistSFNonRigorous(interval A,interval B,int N_search,int max_n);
bool diffusionInNonTwistSF(interval a,interval b,int N_search,int &max_n);
bool diffusionInNonTwistSF(DMap &f,IMap &F,interval a,interval b,int N_search,int max_n);
#endif