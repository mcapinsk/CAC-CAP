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

/////////////////////////////

localMap::localMap(IVector x0_,IVector x1_,IMatrix A_,IMatrix Ainv_,IDissipativeMap &f_)
{
	x0=x0_; x1=x1_;
	A=A_; Ainv=Ainv_;
	f=&f_;
}

/////////////////////////////
// This routine computes a matrix A such that:
// for
//    x = FixedPoint
// and
//    D = Df(x)
// the matrix 
//    A^{-1}*D*A
// is diagonal, and the unstable local coordinate is on the x-axis. 
DMatrix getLocalCoordinates(DDissipativeMap &f,int i)
{
	DVector x=f.FixedPoint(i);
	DMatrix D=f[x];
	DVector rE(2), iE(2);
	DMatrix rVec(2,2), iVec(2,2);
	computeEigenvaluesAndEigenvectors(D,rE,iE,rVec,iVec);
	return rVec;
}

DVector unstableVector(DDissipativeMap &f,int i)
{
	DMatrix A=getLocalCoordinates(f,i);
	DVector v(2);
	v[0]=A[0][0];
	v[1]=A[1][0];
	return v;
}

// Here we take a fixed point 
//    x = FixedPoint
// we consider its unstable eigevector 
//    v=unstableVector
// we take 
//    vl[0] = x-r*v
//    vr[0] = x+r*v
// we compute 
//    vl[i+1] = f(vl[i])
//    vr[i+1] = f(vr[i])
// until we exit the required domain (go below or above B).
vector<DVector> outTrajectory(DDissipativeMap &f,int maxIteratesUpDown,double B,double r,int i,int &side,int UpDown)
{
	DVector x=f.FixedPoint(i);

	vector<DVector> vl,vr;

	DVector v=unstableVector(f,i);

	vl.push_back(x-r*v); 
	vr.push_back(x+r*v);
	
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
vector<DVector> upTrajectory(DDissipativeMap &f,int maxIteratesUpDown,double B,double r,int i,int &side)
{
	return outTrajectory(f,maxIteratesUpDown,B,r,i,side,1);
}
// This computes the guess for a trajectory that will go down
vector<DVector> downTrajectory(DDissipativeMap &f,int maxIteratesUpDown,double B,double r,int i,int &side)
{
	return outTrajectory(f,maxIteratesUpDown,B,r,i,side,-1);
}

////////////////////////////
// Here we check cone conditions.
// The cone C is an initial guess for a cone that will map to itself. 
// The function refines the initial guess for the cone (1,[-L,L]). 
// If the function returns 1, then C will store a sharper cone than the 
// initial guess (1,[-L,L]), for which we have cone conditions.
bool validateCone(localMap &f,interval r,interval L,IVector &C)
{
	IVector v({interval(1),L*interval(-1,1)}); // this is the initial guess for a cone
	IVector x(2);
	// here we iterate  the image of the cone ten times by the derivative
	// to refine the initial guess.
	for(int i=0;i<10;i++) 
	{
		x=r*interval(-1,1)*v;
		v=f[x]*v;
		v[1]=interval(-1,1)*v[1]/v[0].left();
		v[0]=interval(1);
	}

	// We enlarge the guess:
	v[1]=1.1*v[1];
	v[0]=interval(1);
	// We take a neighbourhood which contains the cone with radius r:
	x=r*interval(-1,1)*v;

	// we compute the bound for the image of the cone
	IVector w=f[x]*v;
	// We normalise the result
	w[1]=interval(-1,1)*w[1]/w[0].left();

	// we check cone conditions
	if(! subsetInterior(w[1],v[1])) 
	{
		return 0; // this means that cone conditions failed.
	}
	// we check expansion inside of a cone:
	if(! (abs(w[0])>interval(1.0)))
	{
		return 0;
	}
	C=v; // the C will now pass the validated cone outside of the function
	return 1; // success
}



bool validateChaos(IDissipativeMap &F,DDissipativeMap &Fd,chaosProofParameters par,int &NofIterates,int &failureReason)
{
	interval r=par.r;
	interval B=par.B;
	interval L=par.L;
	double rho=par.rho;
	int maxIteratesUpDown=par.maxIteratesUpDown;

	// While performing tests we wanted to know how many iterates we used.
	// This was computed just out of curiosity. The computed length of the 
	// validated trajectory is not needed for the proof.
	NofIterates=0; 

	// We check if we go from fixedPoint(i) above and below for i=0,1
	for(int i=0;i<2;i++)
	{
		int sideU=0, sideD=0;
		interval h0=rho*r;
		// computation of an initial guess of a trajectory that goes up:
		vector<DVector> qU=upTrajectory(Fd,maxIteratesUpDown,toDouble(B),toDouble(h0),i,sideU);
		if(sideU==0) // this means that there is no trajectory shorter than maxIteratesUpDown that goes up
		{
			if((failureReason==0) or (failureReason==1)) failureReason=1;
			return 0;
		}

		// computation of an initial guess of a trajectory that goes down:
		vector<DVector> qD=downTrajectory(Fd,maxIteratesUpDown,-toDouble(B),toDouble(h0),i,sideD);
		if(sideD==0) // if this is the case we have failed
		{
			if((failureReason==0) or (failureReason==1)) failureReason=1;
			return 0;
		}
		
		// We check the largest number of iterates needed. This is just out of curiosity.
		// This is not essential for the proof.
		if(NofIterates<qU.size()-1) NofIterates=qU.size()-1;
		if(NofIterates<qD.size()-1) NofIterates=qD.size()-1;

		// We compute the local coordinates dor the local map
		IMatrix A=toInterval(getLocalCoordinates(Fd,i));
		IVector x0=F.FixedPoint(i);
		IVector x1=F(x0); 
		IMatrix Ainv=gaussInverseMatrix(A);

		localMap f(x0,x1,A,Ainv,F);
		IVector v0(2);
		if(validateCone(f,r,L,v0)==0)
		{
			failureReason=2;
			return 0;
		}
		// v0 is of the form v0=(v0_x,v0_y)=(1,[-c,c])
		// this means that 
		v0=f.get_A()*v0;

		// We point the cones in the directions of the trajectories:
		IVector vU=sideU*v0;
		IVector vD=sideD*v0;

		IVector q0=F.FixedPoint(i);

		vector<IVector> qU_validated, qD_validated;

		if(validateTrajectoryDown(qD,qD_validated,-B,*(F.map),q0,vD,h0)==0)
		{
			// the "failureReason" was added in order to tweak parameters.
			// Knowing why things failed is not needed for the proof.
			failureReason=3; 
			return 0;
		}
		// we make sure that q0+h0*vD is inside of the cone q0+[0,r]*vD
		if(not(r>h0))
		{
			failureReason=4;
			return 0;
		} 

		h0 = rho*r;
		if(validateTrajectoryUp(qU,qU_validated,B,*(F.map),q0,vU,h0)==0)
		{
			failureReason=3;
			return 0;
		}
		// we make sure that q0+h0*vU is inside of the cone q0+[0,r]*vU
		if(not(r>h0))
		{
			failureReason=4;
			return 0;
		}
	}
	return 1;
}

bool validateChaosStandardMap(IStandardMap &F,DStandardMap &Fd,chaosProofParameters par,int &NofIterates,int &failureReason)
{
	int kappaMax=0.9999*F.get_a().leftBound()/(1.0-F.get_b().leftBound());
	// We try to validate the needed conditions for various kappa:
	for(int kappa=1;kappa<=kappaMax;kappa++)
	{
		F.set_kappa(kappa);
		Fd.set_kappa(kappa);
		if(validateChaos(F,Fd,par,NofIterates,failureReason)==1) return 1;
	}
	return 0;
}

bool validateChaosStandardMap(IStandardMap &F,DStandardMap &Fd,chaosProofParameters par)
{
	int NofIterates, failureReason;
	return validateChaosStandardMap(F,Fd,par,NofIterates,failureReason);
}

// We validate the area in the parameter box a \times b for which we have chaos.
interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Try,int Debth)
{
	if(Try-1==Debth) return interval(0); // this means that we no longer check. 
	F.set_a(a);
	F.set_b(b);
	Fd.set_a(toDouble(a));
	Fd.set_b(toDouble(b));
	par.B=(interval(1.0)/(interval(1.0)-b)).right();

	if(validateChaosStandardMap(F,Fd,par)==1)
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
				Area = Area + chaoticArea(F,Fd,alpha,beta,par,Try+1,Debth);
			} 
		}
		return Area;
	}
}

interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Debth)
{
	// we start with zero attempts, as the chaoticArea() is recursively
	// iterated the Try=0 will be increased until it reaches the max number of iterations
	// indicated by Debth.
	return chaoticArea(F,Fd,a,b,par,0,Debth);
}







