#include "dissipativeMap.h"
#include "validatedTrajectory.h"
#include "utils.h"

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
	setParameter("a",a_);
	a=a_;
}

void DStandardMap::set_b(double b_)
{
	setParameter("b",b_);
	b=b_;
}

DVector DStandardMap::FixedPoint(int i)
{
	double x=0.5*asin(-kappa*(b-1.0)/a)/PI;
	double y=kappa/a;
	if(i==0) return DVector({x,y});
	return DVector({x+0.5,-y});
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
	interval x=0.5*asin(-kappa*(b-1.0)/a)/PI;
	interval y=kappa/a;
	if(i==0) return IVector({x,y});
	return IVector({x+0.5,-y});
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
// until we exit the required domein (go below or above B).
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
	side=0;
	return w;
}

// This computes the guess for a trajectory that will go up/down
vector<DVector> upTrajectory(DDissipativeMap &f,int maxIteratesUpDown,double B,double r,int i,int &side)
{
	return outTrajectory(f,maxIteratesUpDown,B,r,i,side,1);
}
vector<DVector> downTrajectory(DDissipativeMap &f,int maxIteratesUpDown,double B,double r,int i,int &side)
{
	return outTrajectory(f,maxIteratesUpDown,B,r,i,side,-1);
}

////////////////////////////

bool validateCone(localMap &f,interval r,interval L,IVector &C)
{
	IVector v({interval(1),L*interval(-1,1)});
	IVector x(2);
	for(int i=0;i<10;i++)
	{
		x=r*interval(-1,1)*v;
		v=f[x]*v;
		v[1]=interval(-1,1)*v[1]/v[0].left();
		v[0]=interval(1);
	}
	v[1]=1.1*v[1];
	x=r*interval(-1,1)*v;

	IVector w=f[x]*v;
	w[1]=interval(-1,1)*w[1]/w[0].left();

	if(! subsetInterior(w[1],v[1])) 
	{
		return 0;
	}
	// checking expansion inside of a cone:
	if(! (abs(w[0])>interval(1.0)))
	{
		return 0;
	}
	C=v;
	return 1;
}

bool validateChaos(IDissipativeMap &F,DDissipativeMap &Fd,chaosProofParameters par,int &NofIterates,int &failureReason)
{
	interval r=par.r;
	interval B=par.B;
	interval L=par.L;
	double rho=par.rho;
	int maxIteratesUpDown=par.maxIteratesUpDown;

	NofIterates=0;
	for(int i=0;i<2;i++)
	{
		int sideU=0, sideD=0;
		interval h0=rho*r;
		vector<DVector> qU=upTrajectory(Fd,maxIteratesUpDown,toDouble(B),toDouble(h0),i,sideU);
		if(sideU==0) // this means that there is no trajectory shorter than maxIteratesUpDown that goes up
		{
			if((failureReason==0) or (failureReason==1)) failureReason=1;
			return 0;
		}

		vector<DVector> qD=downTrajectory(Fd,maxIteratesUpDown,-toDouble(B),toDouble(h0),i,sideD);
		if(sideD==0)
		{
			if((failureReason==0) or (failureReason==1)) failureReason=1;
			return 0;
		}
		
		if(NofIterates<qU.size()-1) NofIterates=qU.size()-1;
		if(NofIterates<qD.size()-1) NofIterates=qD.size()-1;

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

		IVector vU=sideU*v0;
		IVector vD=sideD*v0;

		IVector y0=F.FixedPoint(i);

		vector<IVector> qU_validated, qD_validated;

		if(validateTrajectoryDown(qD,qD_validated,-B,*(F.map),y0,vD,h0)==0)
		{
			failureReason=3;
			return 0;
		}
		// we make sure that y0+h0*v0 is inside of the cone y0+[0,r]*v0
		if(not(r>h0))
		{
			failureReason=4;
			return 0;
		} 

		h0 = rho*r;
		if(validateTrajectoryUp(qU,qU_validated,B,*(F.map),y0,vU,h0)==0)
		{
			failureReason=3;
			return 0;
		}

		if(not(r>h0))
		{
			failureReason=4;
			return 0;
		}
	}
	return 1;
}

bool validateChaos(DDissipativeMap &Fd,int maxIteratesUpDown,double r,double B,int &NofIterates)
{
	NofIterates=0;
	for(int i=0;i<2;i++)
	{
		int sideU=0, sideD=0;
		vector<DVector> qU=upTrajectory(Fd,maxIteratesUpDown,B,r,i,sideU);
		if(sideU==0) // this means that there is no trajectory shorter than maxIteratesUpDown that goes up
		{
			return 0;
		}

		vector<DVector> qD=downTrajectory(Fd,maxIteratesUpDown,-B,r,i,sideD);
		if(sideD==0)
		{
			return 0;
		}
		if(NofIterates<qU.size()-1) NofIterates=qU.size()-1;
		if(NofIterates<qD.size()-1) NofIterates=qD.size()-1;
	}
	
	return 1;
}

bool validateChaosStandardMap(IStandardMap &F,DStandardMap &Fd,chaosProofParameters par,int &NofIterates,int &failureReason)
{
	int kappaMax=0.9999*F.get_a().leftBound()/(1.0-F.get_b().leftBound());
	for(int kappa=0;kappa<=kappaMax;kappa++)
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

bool chaosInStandardMapNonRigorous(DStandardMap &F,int maxIteratesUpDown,double r,int &NofIterates)
{
	int kappaMax=0.9999*1.0/(1.0-F.get_b());
	double B=1.0/(1.0-F.get_b());
	for(int kappa=0;kappa<=kappaMax;kappa++)
	{
		F.set_kappa(kappa);
		if(validateChaos(F,maxIteratesUpDown,r,B,NofIterates)==1) return 1;
	}
	return 0;
}

interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Try,int Debth)
{
	if(Try-1==Debth) return interval(0);
	F.set_a(a);
	F.set_b(b);
	Fd.set_a(toDouble(a));
	Fd.set_b(toDouble(b));
	par.B=(interval(1.0)/(interval(1.0)-b)).right();

	if(validateChaosStandardMap(F,Fd,par)==1)
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
				Area = Area + chaoticArea(F,Fd,alpha,beta,par,Try+1,Debth);
			} 
		}
		return Area;
	}
}

interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Debth)
{
	return chaoticArea(F,Fd,a,b,par,0,Debth);
}







