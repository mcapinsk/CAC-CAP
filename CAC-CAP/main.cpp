#include <iostream>
#include <chrono>
#include <omp.h>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;
#include "validatedTrajectory.h"
#include "conservativeMap.h"
#include "dissipativeMap.h"
#include "nonTwistSF.h"
#include "plots.h"
#include "utils.h"

interval validateChaosInStandardMapCluster(int Debth)
{
	int nA=700;
	int nB=100;

	interval A(3.0,10.0);
	interval B(0.1,0.8);
	
	int N_of_threads=omp_get_max_threads();
	cout << "Number of threads: " << N_of_threads << endl;

	vector<IStandardMap*> F(N_of_threads);
	vector<DStandardMap*> Fd(N_of_threads);
	vector<chaosProofParameters> par(N_of_threads);
	vector<interval> Area(N_of_threads),a(N_of_threads),b(N_of_threads);

	for(int i=0;i<N_of_threads;i++)
	{
		F[i]=new IStandardMap();
		Fd[i]=new DStandardMap();
		par[i].B = interval(0); // level to go above/below (for below we take -B)
		par[i].r = interval(0.005); // how far from the fixed point we use the cone
		par[i].L = interval(0.01); // initial cone slope
		par[i].rho = 0.5;
		par[i].maxIteratesUpDown=40;
	}
	
	int done=0;
	int i;
	cout << "progress (the program will finish at 1.0): " << endl;
	#pragma omp parallel for private(i)
	for(i=0;i<nA;i++)
	{
		int id=omp_get_thread_num();
		a[id]=part(A,nA,i);
		for(int j=0;j<nB;j++)
		{
			b[id]=part(B,nB,j);
			Area[id] = Area[id] + chaoticArea(*(F[id]),*(Fd[id]),a[id],b[id],par[id],Debth);
		}
		if(i % 7 == 0)
		{
			done++;
			cout << done/100.0 << endl;
		} 
	}
	interval total_Area(0);
	for(int k=0;k<N_of_threads;k++) total_Area = total_Area+ Area[k];
	cout << "total area validated        : " << total_Area << endl;
	cout << "percentage of area validated: " << total_Area/((A.right()-A.left())*(B.right()-B.left())) << endl;
	return total_Area/((A.right()-A.left())*(B.right()-B.left()));
}

void chaosInStandardMapNonRigorous(double a,double b)
{
	DStandardMap F(a,b);
	double r=0.00001;
	int nMax=100;
	int n;
	cout << "result   : " << chaosInStandardMapNonRigorous(F,nMax,r,n) << endl;
	cout << "iterates : " << n << endl;
	plotSystem(F,6);
	//plotSystem(F);
}

void validateDiffusionNoTwistSF()
{
	interval A(0.0,1.0);
	interval B(0.0,1.0);

	int Na=1000;
	int Nb=Na;
	int N_search=2000;
	int max_n=300;
	
	int N_of_threads=omp_get_max_threads();
	cout << "Number of threads: " << N_of_threads << endl;
	vector<ofstream> file(N_of_threads);
	vector<interval> a(N_of_threads), b(N_of_threads);
	vector<IMap*> f(N_of_threads);
	vector<DMap*> fd(N_of_threads);
	for(int i=0;i<N_of_threads;i++)
	{
		file[i].open("plots/NTSF_proof_"+to_string(i)+".txt");
		a[i]=interval(0);
		b[i]=interval(0);
		f[i]=new IMap("par:pi,a,b;var:x,y;fun:x+a*(1-(y-b*sin(2*pi*x))^2),y-b*sin(2*pi*x);");
		f[i]->setParameter("pi",interval::pi());
		fd[i]=new DMap("par:pi,a,b;var:x,y;fun:x+a*(1-(y-b*sin(2*pi*x))^2),y-b*sin(2*pi*x);");
		fd[i]->setParameter("pi",4.0*atan(1.0));
	}

	int i;
	#pragma omp parallel for private(i)
	for(i=0;i<Na;i++)
	{
		int id=omp_get_thread_num();
		a[id]=part(A,Na,i);
		
		for(int j=0;j<Nb;j++)
		{
			b[id]=part(B,Nb,j);
			if(diffusionInNonTwistSF(*(fd[id]),*(f[id]),a[id].mid(),b[id].mid(),N_search,max_n)==1) 
			{
				plot(a[id],b[id],file[id]);
			}
		}
	}
}

void plotDiffusionNoTwistSF()
{
	interval A(0.001,1.0);
	interval B(0.001,1.0);

	int Na=100;
	int Nb=Na;
	int N_search=50;
	int max_n=400;
	ofstream file("plots/NTSF.txt");
	for(int i=0;i<Na;i++)
	{
		interval a=part(A,Na,i);
		for(int j=0;j<Nb;j++)
		{
			interval b=part(B,Nb,j);
			if(diffusionInNonTwistSFNonRigorous(a,b,N_search,max_n)==1) 
				//file << a.mid().leftBound() << " " << b.mid().leftBound() << endl;
				plot(a,b,file);
		}
	}
}



int main()
{
	srand((unsigned)time(NULL));
	rand();
	auto start = chrono::high_resolution_clock::now(); // Start time

	ofstream resultFile("results/StandardMapResult.txt");
	resultFile.precision(15);

	cout.precision(10);
	try
	{	
		validateDiffusionNoTwistSF();
		//plotDiffusionNoTwistSF();
		// return 0;


		// validateChaosInConservativeMaps();
		
		// int accuracy=0;
		// cout << "Please choose accuracy for the Standard Map computation." << endl;
		// cout << "  0 gives a short computation and should validate 0.4 of the area. " << endl;
		// cout << "  1 should result in couple of minutes long computation and should validate 0.67 of the area. " << endl;
		// cout << "  2 should take an hour or two on a desktop computer and should validate 0.88 of the area. " << endl;
		// cout << "  3 should validate 0.98 of the area. This is best done on a cluster, or overnight on a good desktop computer. " << endl;
		// cout << "your choice (0, 1, 2, 3): " ; cin >> accuracy;

		// if((accuracy>=0) and (accuracy<=3))
		// 	resultFile << validateChaosInStandardMapCluster(accuracy) << endl;
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}

  	auto end = chrono::high_resolution_clock::now(); // End time
  	chrono::duration<double> duration = end - start;
    
    cout << "Execution time: " << duration.count() << " seconds" << endl;
    resultFile << "Execution time: " << duration.count() << " seconds" << endl;
  	
  	return 0;
} 
