#include "nonTwistSF.h"
#include "conservativeMap.h"
#include "validatedTrajectory.h"
#include "utils.h"

string NonTwistSFformula()
{
	return "par:pi,a,b;var:x,y;fun:x+a*(1-(y-b*sin(2*pi*x))^2),y-b*sin(2*pi*x);";
}

interval M(interval a,interval b)
{
	return abs(intervalHull(interval(2)*sqrt(interval(1)+interval(3)/a)+abs(b),interval(2)/(a*abs(b)))).right();
}

bool diffusionInNonTwistSF(DMap &f,IMap &F,interval a,interval b,int N_search,int max_n)
{
	f.setParameter("a",toDouble(a));
	f.setParameter("b",toDouble(b));
	interval B=M(a,b)*1.01;

	vector<DVector> q=findShortestTrajectory(f,B.leftBound(),N_search,max_n);

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

bool diffusionInNonTwistSF(interval a,interval b,int N_search,int &max_n)
{
	DMap f(NonTwistSFformula());
	f.setParameter("pi",4.0*atan(1.0));
	f.setParameter("a",toDouble(a));
	f.setParameter("b",toDouble(b));
	interval B=M(a,b)*1.01;

	vector<DVector> q=findShortestTrajectory(f,B.leftBound(),N_search,max_n);

	if(q.size()>=max_n) return 0;
	
	IMap F(NonTwistSFformula());
	F.setParameter("pi",interval::pi());
	F.setParameter("a",a);
	F.setParameter("b",b);

	IVector q0(2);
	q0[0]=q[0][0];
	q0[1]=-B;
	IVector v0({interval(1),interval(0)});
	
	interval h0;
	vector<IVector> q_validated;
	max_n=q.size();
	return validateTrajectoryUp(q,q_validated,B,F,q0,v0,h0);
}

bool diffusionInNonTwistSFNonRigorous(interval a,interval b,int N_search,int max_n)
{
	DMap f(NonTwistSFformula());
	f.setParameter("pi",4.0*atan(1.0));
	f.setParameter("a",toDouble(a));
	f.setParameter("b",toDouble(b));
	
	interval B=M(a,b)*1.01;

	vector<DVector> q=findShortestTrajectory(f,B.leftBound(),N_search,max_n);
	//cout << q.size() << " " << B.leftBound() << " " << q[0] << q[q.size()-1]<< endl;
	if(q.size()<max_n) return 1;
	return 0;
}

interval diffusionAreaNonTwistSF(interval a,interval b,int N_search,int max_n,int Try,int Debth)
{
	if(Try-1==Debth) return interval(0);

	if(diffusionInNonTwistSF(a,b,N_search,max_n)==1)
	{
		return (a.right()-a.left())*(b.right()-b.left());
	}else
	{
		interval Area(0);
		int N=5;
		for(int i=0;i<N;i++)
		{
			interval alpha=part(a,N,i);
			for(int k=0;k<N;k++)
			{
				interval beta=part(b,N,k);
				Area = Area + diffusionAreaNonTwistSF(alpha,beta,N_search,max_n,Try+1,Debth);
			} 
		}
		return Area;
	}
}