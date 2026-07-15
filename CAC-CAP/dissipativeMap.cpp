#include "dissipativeMap.h"
#include "validatedTrajectory.h"
#include "utils.h"

// The constructor simple sets up the map formula and parameters.
// (During the course of the proof the parameter kappa can be modified
// using set_kappa().)
DStandardMap::DStandardMap(double a_,double b_) : DDissipativeMap("par:pi,a,b;var:x,y;fun:x-a*(y-sin(2*pi*x))/b,(y-sin(2*pi*x))/b;")
{
	a=a_;
	b=b_;
	PI=4.0*atan(1.0);
	kappa=0;
	setParameter("pi",PI);
	setParameter("a",a);
	setParameter("b",b);
}

// we can also initiate the map without specifying the parameters, in which case 
// a and b are set to zero. 
// (During the course of the proof the parameter a and b can be modified by
// using set_a() and set_b() functions.)
DStandardMap::DStandardMap() : DDissipativeMap("par:pi,a,b;var:x,y;fun:x-a*(y-sin(2*pi*x))/b,(y-sin(2*pi*x))/b;")
{
	a=0.0;
	b=0.0;
	PI=4.0*atan(1.0);
	kappa=0;
	setParameter("pi",PI);
	setParameter("a",a);
	setParameter("b",b);
}

void DStandardMap::set_a(double a_)
{
	// this is needed so that the parameter "a" of the map in DDissipativeMap 
	// is updated:
	setParameter("a",a_); 
	a=a_;
}

void DStandardMap::set_b(double b_)
{
	// this is needed so that the parameter "b" of the map in DDissipativeMap 
	// is updated:
	setParameter("b",b_);
	b=b_;
}

DVector DStandardMap::FixedPoint(int i)
{

	double x=0.5*asin(-kappa*(b-1.0)/a)/PI;
	double y=kappa/a;
	if(i==0) return DVector({x,y}); // when i=0 we return (x_{-kappa}, y_{-kappa})
	
	x=0.5*asin(kappa*(b-1.0)/a)/PI;
	y=-kappa/a;
	return DVector({x,y}); // when i=1 we return (x_{kappa}, y_{kappa})
}

////////////////////////////////////////////////
IStandardMap::IStandardMap(interval a_,interval b_) : IDissipativeMap("par:pi,a,b;var:x,y;fun:x-a*(y-sin(2*pi*x))/b,(y-sin(2*pi*x))/b;")
{
	a=a_;
	b=b_;
	PI=interval::pi();
	kappa=0;
	setParameter("pi",PI);
	setParameter("a",a);
	setParameter("b",b);
}

IStandardMap::IStandardMap() : IDissipativeMap("par:pi,a,b;var:x,y;fun:x-a*(y-sin(2*pi*x))/b,(y-sin(2*pi*x))/b;")
{
	a=interval(0);
	b=interval(0);
	PI=interval::pi();
	kappa=0;
	setParameter("pi",PI);
	setParameter("a",a);
	setParameter("b",b);
}

void IStandardMap::set_a(interval a_)
{
	setParameter("a",a_);
	a=a_;
}

void IStandardMap::set_b(interval b_)
{
	setParameter("b",b_);
	b=b_;
}

IVector IStandardMap::FixedPoint(int i)
{
	interval x=interval(0.5)*asin(-kappa*(b-interval(1.0))/a)/PI;
	interval y=kappa/a;
	if(i==0) return IVector({x,y}); // when i=0 we return (x_{-kappa}, y_{-kappa})
	
	x=interval(0.5)*asin(kappa*(b-interval(1.0))/a)/PI;
	y=-kappa/a;
	return IVector({x,y}); // when i=1 we return (x_{kappa}, y_{kappa})
}

DVector unstableVector(DStandardMap &f,int i)
{
	DVector x=f.FixedPoint(i);
	DMatrix D=f[x];
	DVector rE(2), iE(2);
	DMatrix rVec(2,2), iVec(2,2);
	computeEigenvaluesAndEigenvectors(D,rE,iE,rVec,iVec);

	DVector v(2);
	v[0]=rVec[0][0];
	v[1]=rVec[1][0];

	double lambda = rE[0];
	double V=sqrt(v[0]*v[0]+v[1]*v[1]);
	double scaling=0.1/(V*lambda);

	return scaling*v;
}

////////////////////////////////
// Here we consider a sector
//    S = {x0+s*v : s \in [0,h]} 
// and validate that
//    F(S) \subset F(x0) + [-0.5,0.5] x [-1,1]
bool sectorInclusion(IStandardMap &F,IVector v,interval h,int i)
{
	IVector x0=F.FixedPoint(i);

	int kappa = F.get_kappa();
	IVector Fx0(2);

	// if i=0 then x0 = (x_{-kappa}, y_{-kappa})
	// and F(x0) = x0 - (kappa,0)
	if(i==0) Fx0 = x0 - IVector({interval(kappa),0});

	// if i=1 when x0 = (x_{kappa}, y_{kappa})
	// and F(x0) = x0 + (kappa,0)
	if(i==1) Fx0 = x0 + IVector({interval(kappa),0});

	// We choose a set Bkappa
	// which we are sure will be contained in
	// F(x0) + [-0.5,0.5] x [-1,1]
	IVector Bkappa(2);
	Bkappa[0] = intervalHull(Fx0[0].right()-0.49,Fx0[0].left()+0.49);
	Bkappa[1] = intervalHull(Fx0[1].right()-0.99,Fx0[1].left()+0.99);

	// We consider S = {x0+s*v : s \in [0,h]}.
	// We use the mean value theorem:
	//    F(x0+s*v) \in F(x0) + [DF]*v*[0,h]
	// where [DF] is the interval arithmetic bound on DF(S).
	interval H=intervalHull(interval(0),h);
	IVector FS = Fx0 + H*(F[x0+H*v]*v);	

	if(subsetInterior(FS,Bkappa)) return true;
	return false;
}

// Here we take a fixed point 
//    x = FixedPoint
// we consider its unstable eigevector 
//    v=unstableVector
// we take 
//    vl[0] = x-h*v
//    vr[0] = x+h*v
// we compute 
//    vl[i+1] = f(vl[i])
//    vr[i+1] = f(vr[i])
// until we exit the required domain (go below or above B).
vector<DVector> outTrajectory(DStandardMap &f,int maxIteratesUpDown,double B,int i,int &side,int UpDown,double h)
{
	DVector x=f.FixedPoint(i);

	vector<DVector> vl,vr;

	DVector v=unstableVector(f,i);

	vl.push_back(x-h*v); 
	vr.push_back(x+h*v);
	
	for(int i=0;i<maxIteratesUpDown;i++)
	{
		if(UpDown==1)
		{
			if(vl[i][1]>B)
			{ 
				side=-1;
				return vl;
			}
			if(vr[i][1]>B)
			{
				side=1;
				return vr;
			}
		}
		if(UpDown==-1)
		{
			if(vl[i][1]<B)
			{ 
				side=-1;
				return vl;
			}
			if(vr[i][1]<B)
			{
				side=1;
				return vr;
			}
		}
		
		vl.push_back(f(vl[i])); 
		vr.push_back(f(vr[i]));
	}
	vector<DVector> w;
	// by setting side=0 we mark that we have failed to find a candidate for a trajectory 
	// which goes above/below the required level.
	side=0; 
	return w; // we return an empty sequence of vectors.
}

// This computes the guess for a trajectory that will go up
vector<DVector> upTrajectory(DStandardMap &f,int maxIteratesUpDown,double B,int i,int &side,double h)
{
	return outTrajectory(f,maxIteratesUpDown,B,i,side,1,h);
}

// This computes the guess for a trajectory that will go down
vector<DVector> downTrajectory(DStandardMap &f,int maxIteratesUpDown,double B,int i,int &side,double h)
{
	return outTrajectory(f,maxIteratesUpDown,B,i,side,-1,h);
}

// This function validates the conditions (I-VIII) for the 
// StandardMap F. The map F can have parameters (a,b) chosen as intervals.
bool validateChaos(IStandardMap &F,DStandardMap &Fd,chaosProofParameters par,int &NofIterates)
{
	interval B=par.B;
	int maxIteratesUpDown=par.maxIteratesUpDown;

	// We check if we go from fixedPoint(i) above and below for i=0,1
	for(int i=0;i<2;i++)
	{
		interval h=2.0;
		bool success=false;
		for(int j=0;j<3;j++)
		{
			if(success==true) continue;
			
			h=h/2;
			int sideU=0, sideD=0;
			// computation of an initial guess of a trajectory that goes up:
			vector<DVector> qU=upTrajectory(Fd,maxIteratesUpDown,toDouble(B),i,sideU,h.leftBound());
			if(sideU==0) continue; // this means that there is no trajectory shorter than maxIteratesUpDown that goes up
			
			// computation of an initial guess of a trajectory that goes down:
			vector<DVector> qD=downTrajectory(Fd,maxIteratesUpDown,-toDouble(B),i,sideD,h.leftBound());
			if(sideD==0) continue; // if this is the case we have failed
		
			IVector q0=F.FixedPoint(i);
			IVector vU=IVector({qU[0][0],qU[0][1]})-q0;
			IVector vD=IVector({qD[0][0],qD[0][1]})-q0;

			vector<IVector> qU_validated, qD_validated;

			interval h0=1.0;

			// Below, validateTrajectoryDown() computes h0 and validates that a trajectory starting from
			//    q0+h0*vD
			// goes below -B.
			if(validateTrajectoryDown(qD,qD_validated,-B,*(F.map),q0,vD,h0)==0)	continue;
			// We consider
			//    S = {q0+s*vD : s \in [0,h0]}
			// and validate that F(S) \subset F(q0) + [-0.5,0.5] x [-1,1].
			if(sectorInclusion(F,vD,h0,i)==0) continue;

			h0 = 1.0;

			// Below, validateTrajectoryUp() computes h0 and validates that a trajectory starting from
			//    q0+h0*vU
			// goes above B.
			if(validateTrajectoryUp(qU,qU_validated,B,*(F.map),q0,vU,h0)==0) continue;
			// We consider
			//    S = {q0+s*vU : s \in [0,h0]}
			// and validate that F(S) \subset F(q0) + [-0.5,0.5] x [-1,1].
			if(sectorInclusion(F,vU,h0,i)==0) continue;

			// We check the largest number of iterates needed. This is just out of curiosity.
			// This is not essential for the proof.
			if(NofIterates<(int)qU.size()-1) NofIterates=qU.size()-1;
			if(NofIterates<(int)qD.size()-1) NofIterates=qD.size()-1;
			success = true;
		}
		if(success==false) return 0;
	}

	return 1;
}

bool validateChaosStandardMap(IStandardMap &F,DStandardMap &Fd,chaosProofParameters par,int &NofIterates)
{
	int kappaMax=0.9999*F.get_a().leftBound()/(1.0-F.get_b().leftBound());
	// We try to validate the needed conditions for various kappa:
	for(int kappa=1;kappa<=kappaMax;kappa++)
	{
		F.set_kappa(kappa);
		Fd.set_kappa(kappa);
		if(validateChaos(F,Fd,par,NofIterates)==1) return 1;
	}
	return 0;
}

// We validate the area in the parameter box a \times b for which we have chaos.
interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Try,int Debth,int &NofIterates)
{
	if(Try-1==Debth) return interval(0); // this means that we no longer check. 
	F.set_a(a);
	F.set_b(b);
	Fd.set_a(toDouble(a));
	Fd.set_b(toDouble(b));
	par.B=(interval(1.0)/(interval(1.0)-b)).right();

	if(validateChaosStandardMap(F,Fd,par,NofIterates)==1)
	{
		return (a.right()-a.left())*(b.right()-b.left());
	}else // we subdivide the box a \times b into 25 pieces and try again on each of them:
	{
		interval Area(0);
		int N=5;
		for(int i=0;i<N;i++)
		{
			interval alpha=part(a,N,i);
			for(int k=0;k<N;k++)
			{
				interval beta=part(b,N,k);
				Area = Area + chaoticArea(F,Fd,alpha,beta,par,Try+1,Debth,NofIterates);
			} 
		}
		return Area;
	}
}

interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Debth,int &NofIterates)
{
	// we start with zero attempts, as the chaoticArea() is recursively
	// iterated the Try=0 will be increased until it reaches the max number of iterations
	// indicated by Debth.
	return chaoticArea(F,Fd,a,b,par,0,Debth,NofIterates);
}







