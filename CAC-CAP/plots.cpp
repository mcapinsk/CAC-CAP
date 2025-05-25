#include "plots.h"
#include "utils.h"

// Selecting only orbits which stay within y\in [-B,B]
// for M interates.
void plotSystem(DStandardMap &F,int M)
{	
	int k=0;
	while(k<30)
	{
		double U1=uniformSapmle();
		double U2=uniformSapmle();

		double B=1.0/(1.0-F.get_b());
		
		DVector x0({U1,2*B*U2-B});
		DVector x=x0;
		ofstream file("plots/sym"+to_string(k)+".txt");
		int N=0;
		while((abs(x[1])<3*B) and (N<M+1)) 
		{
			x=F(x);
			N++;
		}
		if(N>M)
		{
			x=x0;
			for(int i=0;i<M+1;i++)
			{
				file << x[0] << " " << x[1] << endl;
				x=F(x);
			}
			file << x[0] << " " << x[1] << endl;
			k++;
		}
	}
}

// selecting random trajectories
void plotSystem(DStandardMap &F)
{	
	int k=0;
	while(k<30)
	{
		double U1=uniformSapmle();
		double U2=uniformSapmle();

		double B=1.0/(1.0-F.get_b());
		B=1.0;
		DVector x0({U1,2*B*U2-B});
		DVector x=x0;
		ofstream file("plots/sym"+to_string(k)+".txt");
		int N=0;
		while(abs(x[1])<3*B) 
		{
			x=F(x);
			N++;
		}
		x=x0;
		for(int i=0;i<N;i++)
		{
			file << x[0] << " " << x[1] << endl;
			x=F(x);
		}
		file << x[0] << " " << x[1] << endl;
		k++;
	}
}