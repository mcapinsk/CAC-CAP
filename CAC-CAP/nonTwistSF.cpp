#include "nonTwistSF.h"
#include "conservativeMap.h"
#include "validatedTrajectory.h"
#include "utils.h"

interval M(interval a,interval b)
{
	// max{2*sqrt(1+3/a)+|b|,2/(a*|b|)}
	return abs(intervalHull(interval(2)*sqrt(interval(1)+interval(3)/a)+abs(b),interval(2)/(a*abs(b)))).right();
}

bool diffusionInNonTwistSF(DMap &f,IMap &F,interval a,interval b,int N_search,int max_n)
{
	f.setParameter("a",toDouble(a));
	f.setParameter("b",toDouble(b));
	interval B=M(a,b)*1.01;

	// We find a non-rogorous candidate for a trajectory which goes up
	vector<DVector> q=findShortestTrajectoryUp(f,B.leftBound(),N_search,max_n);

	if(q.size()>=max_n) return 0;
	
	F.setParameter("a",a);
	F.setParameter("b",b);

	IVector q0(2);
	q0[0]=q[0][0];
	q0[1]=-B;
	IVector v0({interval(1),interval(0)});
	
	interval h0;
	vector<IVector> q_validated;
	return validateTrajectoryUp(q,q_validated,B,F,q0,v0,h0);
}
