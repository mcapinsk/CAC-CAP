#include "validatedTrajectory.h"
using namespace capd::matrixAlgorithms;

IVector combine(vector<IVector> x)
{
	// last vector x[n]={h0,h1}
	int n=x.size()-1;
	IVector X(2*n+2);
	for(int i=0;i<=n;i++)
	{
		X[2*i]  =x[i][0];
		X[2*i+1]=x[i][1];
	}
	return X;
}
vector<IVector> split(IVector X)
{
	// last vector x[n]={h0,h1}
	int n=X.dimension()/2-1;
	vector<IVector> x(n+1);
	for(int i=0;i<=n;i++)
	{
		x[i]=IVector(2);
		x[i][0]=X[2*i];
		x[i][1]=X[2*i+1];
	}
	return x;
}

// This function computes:
//    r[0] = f(x[0])-x[1]
//    r[1] = f(x[1])-x[2]
//     ...
//    r[n-2] = f(x[n-2])-x[n-1]
//    r[n-1] = f(x[n-1])-(y[1]+h1*v[1])
//    r[n]   = f(y[0]+h0*v[0])-x[0]
// and returns r.
IVector F(IVector X,IMap &f,vector<IVector> y,vector<IVector> v)
{
	vector<IVector> x=split(X);
	int n=x.size()-1;
	interval h0 = x[n][0];
	interval h1 = x[n][1];
	vector<IVector> r(n+1);

	for(int i=0;i<n-1;i++)
	{
		r[i] = f(x[i])-x[i+1];
	}
	r[n-1]=f(x[n-1])-(y[1]+h1*v[1]);
	r[n]=f(y[0]+h0*v[0])-x[0];
	return combine(r);
}

IMatrix DF(IVector X,IMap &f,vector<IVector> y,vector<IVector> v)
{
	vector<IVector> x=split(X);
	int n=x.size()-1;
	interval h0 = x[n][0];
	interval h1 = x[n][1];
	IMatrix dF(2*n+2,2*n+2);
	for(int i=0;i<n;i++)
	{
		IMatrix D=f[x[i]];
		dF[2*i  ][2*i] = D[0][0]; dF[2*i  ][2*i+1] = D[0][1];
		dF[2*i+1][2*i] = D[1][0]; dF[2*i+1][2*i+1] = D[1][1];
	}
	for(int i=0;i<2*n-2;i++) dF[i][i+2]=-1;

	dF[2*n-2][2*n+1]=-v[1][0];
	dF[2*n-1][2*n+1]=-v[1][1];

	dF[2*n  ][0]=-1; dF[2*n  ][1]=0;
	dF[2*n+1][0]= 0; dF[2*n+1][1]=-1;

	IVector w=f[y[0]+h0*v[0]]*v[0];
	dF[2*n  ][2*n]=w[0];
	dF[2*n+1][2*n]=w[1];
	return dF;
}

IMatrix getInverseGuess(IMatrix A)
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
//    q[0] = y[0]+h0*v[0]
//    q[1] = f(q[0])
//    q[2] = f(q[1])
//    q[3] = f(q[2])
//     ...
//    q[n] = f(q[n-1])
//	  q[n+1] = f(q[n]) = y[1]+h1*v[1]
// which is contained in a sequence of boxes
//    q[0]   \in p[0] 
//       ...
//    q[n+1] \in p[n+1]
// The routine computes:
//    vector<IVector> &p
//    interval &h0
// and returns 1 if they are computed correctly.
bool validatedTrajectory(vector<IVector> &p,IVector y0,vector<IVector> v,interval &h0,IMap &f)
{
	int n=p.size()-2;
	vector<IVector> y(2);
	y[0]=y0;
	y[1]=p[n+1];

	vector<IVector> x(n+1);
	for(int i=0;i<n;i++) x[i]=p[i+1];
	x[n]=IVector(2); // this is {h0,h1}
	x[n][0]=h0; // this is the initial guess for h0. It will be validated and potentially improved by the routine. 
	x[n][1]=0.0; // this will be computed by the routine

	IVector X=combine(x);
	IVector X0=midVector(X);

	// we compute the matrix C needed for the Krawczyk operator. 
	// The C will remain fixed.
	IMatrix C=getInverseGuess(DF(X0,f,y,v));
	
	
	int m=C.numberOfColumns();
	IMatrix Id(m,m);
	for(int i=0;i<m;i++) Id[i][i]=interval(1.0);

	// blowup of the point X:
	for(int i=0;i<3;i++)
	{
		X=X0 - C*F(X0,f,y,v) + (Id - C*DF(X,f,y,v))*(X-X0);
		X0=midVector(X);
	}
	X=X0+2.0*(X-X0);

	// this is the heart of the code:
	// the Krawczyk method:
	IVector K = X0 - C*F(X0,f,y,v) + (Id - C*DF(X,f,y,v))*(X-X0);

	if(subsetInterior(K,X)) // :)
	{
		x=split(X);
		h0=x[n][0]; // this modifies the h0 passed to the function
		interval h1=x[n][1];
		p[0] = y[0] + h0*v[0];
		p[n+1] = y[1] + h1*v[1];
		for(int i=0;i<n;i++) p[i+1]=x[i];
		return 1;
	}
	// cout << "X size:" << X-midVector(X) << endl;
	// cout << "K size:" << K-midVector(K) << endl;
	// cout << "X - K :" << X-K << endl;
	
	return 0;
}

// This routine does the same thing as the one above,
// the only difference is that we compute the vector v1 by trial and
// error. We check if v[1] = (0,1) works, if yes great, if not we try 
// v[1] = (1,0).
bool validatedTrajectory(vector<IVector> &p,IVector y0,IVector v0, interval &h0,IMap &f)
{
	vector<IVector> v(2);
	v[0] = v0;
	
	v[1] = IVector({1,0});
	vector<IVector> q=p;

	// Remark:
	// by calling validatedTrajectory(q,y0,v,h0,f)
	// the routine computes q and h0.
	if(validatedTrajectory(q,y0,v,h0,f)==0)
	{
		v[1] = IVector({0,1});
		return validatedTrajectory(p,y0,v,h0,f);
	}
	p=q;
	return 1;
}

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

bool validateTrajectoryUp(vector<DVector> q_guess,vector<IVector> &q_validated,interval L,IMap &F,IVector y0,IVector v0,interval &h0)
{
	int n=q_guess.size();
	q_validated=to_interval(q_guess);
	bool res=validatedTrajectory(q_validated,y0,v0,h0,F);
	if(res==0) return res;
	if(not (q_validated[n-1][1]>L)) return 0;
	return 1;
}

bool validateTrajectoryDown(vector<DVector> q_guess,vector<IVector> &q_validated,interval L,IMap &F,IVector y0,IVector v0,interval &h0)
{
	int n=q_guess.size();
	q_validated=to_interval(q_guess);
	bool res=validatedTrajectory(q_validated,y0,v0,h0,F);
	if(res==0) return res;
	if(not (q_validated[n-1][1]<L)) return 0;
	return 1;
}

