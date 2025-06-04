#include "validatedTrajectory.h"
using namespace capd::matrixAlgorithms;

// This routine takes a sequence of vercors n in R^2 and
// stacks them up in a vector X in R^(2*n).
IVector combine(vector<IVector> x)
{
	int n=x.size();
	IVector X(2*n);
	for(int i=0;i<n;i++)
	{
		X[2*i]  =x[i][0];
		X[2*i+1]=x[i][1];
	}
	return X;
}

// This routine has the opposite role to the one above.
vector<IVector> split(IVector X)
{
	int n=X.dimension()/2;
	vector<IVector> x(n);
	for(int i=0;i<n;i++)
	{
		x[i]=IVector(2);
		x[i][0]=X[2*i];
		x[i][1]=X[2*i+1];
	}
	return x;
}

/////////////////////////////////////////////
// This function computes:
//    r[0] = f(p[0])-p[1]
//    r[1] = f(p[1])-p[2]
//     ...
//    r[n-2] = f(p[n-2])-p[n-1]
//    r[n-1] = f(p[n-1])-(q[1]+h1*v[1])
//    r[n]   = f(q[0]+h0*v[0])-p[0]
// and returns r. See (8) from the paper.
//
// The convention is that the variables h0 and h1 are the last two
// coefficients of the vector X.
IVector F(IVector X,IMap &f,vector<IVector> q,vector<IVector> v)
{
	vector<IVector> p=split(X);
	int n=p.size()-1;
	interval h0 = p[n][0];
	interval h1 = p[n][1];
	vector<IVector> r(n+1);

	for(int i=0;i<n-1;i++)
	{
		r[i] = f(p[i])-p[i+1];
	}
	r[n-1]=f(p[n-1])-(q[1]+h1*v[1]);
	r[n]=f(q[0]+h0*v[0])-p[0];
	return combine(r);
}

// This is the derivative of the above map F.

IMatrix DF(IVector X,IMap &f,vector<IVector> q,vector<IVector> v)
{
	vector<IVector> p=split(X);
	int n=p.size()-1;
	interval h0 = p[n][0];
	interval h1 = p[n][1];
	IMatrix dF(2*n+2,2*n+2);
	for(int i=0;i<n;i++)
	{
		IMatrix D=f[p[i]];
		dF[2*i  ][2*i] = D[0][0]; dF[2*i  ][2*i+1] = D[0][1];
		dF[2*i+1][2*i] = D[1][0]; dF[2*i+1][2*i+1] = D[1][1];
	}
	for(int i=0;i<2*n-2;i++) dF[i][i+2]=-1;

	dF[2*n-2][2*n+1]=-v[1][0];
	dF[2*n-1][2*n+1]=-v[1][1];

	dF[2*n  ][0]=-1; dF[2*n  ][1]=0;
	dF[2*n+1][0]= 0; dF[2*n+1][1]=-1;

	IVector w=f[q[0]+h0*v[0]]*v[0];
	dF[2*n  ][2*n]=w[0];
	dF[2*n+1][2*n]=w[1];
	return dF;
}

// This function computes a guess of an inverse of a matrix A.
// This is not an interval arithmetic computation. This approximate
// inverse is not a rigorous enclosure of the inverse of A. It does
// not need to be. This approximation is later used as the matrix
// used for the Krawczyk method. Krawczyk method can work with any
// matrix, as long as the assumptions of the theorem can be validated.
IMatrix getCloseToInverseGuess(IMatrix A)
{
	int n=A.numberOfColumns();
	DMatrix C(n,n);
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++) C[i][k]=A[i][k].mid().leftBound();
	}
	C=gaussInverseMatrix(C);

	for(int i=0;i<n;i++)
	{
		for(int k=0;k<n;k++) A[i][k]=C[i][k];
	}
	return A;
}

// This routine validates that there exists a trajectory 
//    p[0] = q[0]+h0*v[0]
//    p[1] = f(p[0])
//    p[2] = f(p[1])
//    p[3] = f(p[2])
//     ...
//    p[n] = f(p[n-1])
//	  p[n+1] = f(p[n]) = q[1]+h1*v[1]
// which is contained in a sequence of boxes
//    p[0]   \in P[0] 
//       ...
//    p[n+1] \in P[n+1]
// The routine computes:
//    vector<IVector> &P
//    interval &h0
// and returns 1 if they are computed correctly.
//
// The vector<IVector> P is given on input and it should consist of an initial guess
// on the trajectory leading from q0+h0*v[0]. The last coefficient of P should consist
// of the final point of the trajectory q1. This means that the final point 
// after validating will be q1+h1*v[1], for some h1 close to zero. 
//
// If the validation is succesful, the validated enclosure of the trajectory is
// stored in P. The function also modifies h0 to ensure that q0+h0*v[0] is 
// a validated enclosure of the starting point.
bool validatedTrajectory(vector<IVector> &P,IVector q0,vector<IVector> v,interval &h0,IMap &f)
{
	int n=P.size()-2;
	vector<IVector> q(2);
	q[0]=q0;
	q[1]=P[n+1];

	vector<IVector> p(n+1);
	for(int i=0;i<n;i++) p[i]=P[i+1];
	p[n]=IVector(2); // this is {h0,h1}
	p[n][0]=h0; // this is the initial guess for h0. It will be validated and potentially improved by the routine. 
	p[n][1]=0.0; // this will be computed by the routine

	IVector X=combine(p);
	IVector X0=midVector(X);

	// we compute the matrix C needed for the Krawczyk operator. 
	// The C will remain fixed.
	IMatrix C=getCloseToInverseGuess(DF(X0,f,q,v));
	
	
	int m=C.numberOfColumns();
	IMatrix Id(m,m);
	for(int i=0;i<m;i++) Id[i][i]=interval(1.0);

	// blowup of the point X:
	for(int i=0;i<3;i++)
	{
		X=X0 - C*F(X0,f,q,v) + (Id - C*DF(X,f,q,v))*(X-X0);
		X0=midVector(X);
	}
	X=X0+2.0*(X-X0);

	// this is the heart of the code:
	// the Krawczyk method:
	IVector K = X0 - C*F(X0,f,q,v) + (Id - C*DF(X,f,q,v))*(X-X0);

	if(subsetInterior(K,X)) // If this is true, the validation is successful.
	{
		p=split(X);
		h0=p[n][0]; // this modifies the h0 passed to the function
		interval h1=p[n][1];
		P[0] = q[0] + h0*v[0]; // this is the starting point of the trajectory
		P[n+1] = q[1] + h1*v[1]; // this is the final point of the trajectory
		for(int i=0;i<n;i++) P[i+1]=p[i];
		return 1;
	}
	
	return 0;
}

// This routine does the same thing as the one above,
// the only difference is that we compute the vector v1 by trial and
// error. We check if v[1] = (0,1) works, if yes great, if not we try 
// v[1] = (1,0).
bool validatedTrajectory(vector<IVector> &P,IVector y0,IVector v0, interval &h0,IMap &f)
{
	vector<IVector> v(2);
	v[0] = v0;
	
	v[1] = IVector({1,0});
	vector<IVector> q=P;

	// Remark:
	// by calling validatedTrajectory(q,y0,v,h0,f)
	// the routine computes q and h0.
	if(validatedTrajectory(q,y0,v,h0,f)==0)
	{
		v[1] = IVector({0,1});
		return validatedTrajectory(P,y0,v,h0,f);
	}
	P=q;
	return 1;
}

// This is a technical function to convert a DVector to IVector.
// Its only purpose is to convert a sequence of DVectors which 
// stores a non-rigorous guess of the trajectory to a sequence of IVectors,
// which are then passed to the validatedTrajectory() function.
vector<IVector> to_interval(vector<DVector> p)
{
	int n=p.size();
	vector<IVector> q(n);
	for(int i=0;i<n;i++)
	{
		q[i]=IVector({p[i][0],p[i][1]});
	}
	return q;
}

// This function validates that a trajectory which starts from
// q0+h0*v0
// goes above L.
// 
// The function stores the validated trajectory in p_validated, and also modifies
// h0 so that 
// q0+h0*v0
// is a true enclosure of a starting point which goes above L.
bool validateTrajectoryUp(vector<DVector> p_guess,vector<IVector> &p_validated,interval L,IMap &F,IVector q0,IVector v0,interval &h0)
{
	int n=p_guess.size();
	p_validated=to_interval(p_guess);
	bool res=validatedTrajectory(p_validated,q0,v0,h0,F);
	if(res==0) return res;
	if(not (p_validated[n-1][1]>L)) return 0;
	return 1;
}

// This function performs the same task as the one above,
// but checks if the trajectory goes below L.
bool validateTrajectoryDown(vector<DVector> p_guess,vector<IVector> &p_validated,interval L,IMap &F,IVector q0,IVector v0,interval &h0)
{
	int n=p_guess.size();
	p_validated=to_interval(p_guess);
	bool res=validatedTrajectory(p_validated,q0,v0,h0,F);
	if(res==0) return res;
	if(not (p_validated[n-1][1]<L)) return 0;
	return 1;
}

