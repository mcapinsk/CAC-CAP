#ifndef validatedTrajectory_h
#define validatedTrajectory_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;

//////////////////////////////
// This routine validates a rigorous enclosure of a trajectory with initial point:
//   q0+h0*v0
// which goes above the level L.
// The validation is performed by using the Krawczyk method described in the paper.
// The q_guess is a non-rigorously computed guess for where the trajectory is.
// The routine validates that close to the q_guess we have the true trajectory,
// which is then written into the sequence of interval vectors "p_validated".
// The function also modifies h0 so that 
//   q0+h0*v0
// is a true enclosure of a starting point which goes above L.
bool validateTrajectoryUp(vector<DVector> p_guess,vector<IVector> &p_validated,interval L,IMap &F,IVector q0,IVector v0,interval &h0);

//////////////////////////////
// This routine validates a rigorous enclosure of a trajectory with initial point:
//   q0+h0*v0
// which goes below the level L.
bool validateTrajectoryDown(vector<DVector> p_guess,vector<IVector> &p_validated,interval L,IMap &F,IVector q0,IVector v0,interval &h0);


#endif