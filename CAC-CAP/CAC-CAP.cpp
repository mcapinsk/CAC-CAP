#include <iostream>
#include <chrono>
#include <omp.h>
#include <cstdio>
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
#include "utils.h"

void producePlotFile(int n_of_threads)
{
	ofstream file("thm_1_3/plotNTSF.txt");
	file << "set xrange [0:1]" << endl;
	file << "set xtics 0,0.2,1" << endl;
	file << "set yrange [0:1]" << endl;
	file << "set ytics 0,0.2,1 " << endl << endl;

	file << "plot \\" << endl;
	for(int i=0;i<n_of_threads;i++)
	{
		file << "'NTSF_proof_"+to_string(i)+".txt' lw 1.3 lc rgb 'red' w d t '', \\" << endl;
	}
	file << endl << "set terminal postscript eps enhanced size 2.0in,2.0in color font 'Helvetica,12' " << endl;
	file << "set output 'Fig-ntsf.eps'" << endl;
	file << "replot" << endl;
	file << "set terminal qt" << endl;
	file << "reset" << endl;

	// removing aur (authors) results which are downloaded from github:
	for(int i=n_of_threads;i<48;i++)
	{
		string filename="thm_1_3/NTSF_proof_"+to_string(i)+".txt";
		std::remove(filename.c_str());
	}
}

void validateDiffusionNoTwistSF(int Mesh_Size)
{
	interval A(0.0,1.0);
	interval B(0.0,1.0);

	int Na=Mesh_Size;
	int Nb=Na;
	int N_search=2000;
	int max_n=300;
	
	int N_of_threads=omp_get_max_threads();
	cout << "Number of threads: " << N_of_threads << endl;
	producePlotFile(N_of_threads);
	vector<ofstream> file(N_of_threads);
	vector<interval> a(N_of_threads), b(N_of_threads);
	vector<IMap*> f(N_of_threads);
	vector<DMap*> fd(N_of_threads);
	vector<int> counter(N_of_threads);
	for(int i=0;i<N_of_threads;i++)
	{
		file[i].open("thm_1_3/NTSF_proof_"+to_string(i)+".txt");
		a[i]=interval(0);
		b[i]=interval(0);
		f[i]=new IMap("par:pi,a,b;var:x,y;fun:x+a*(1-(y-b*sin(2*pi*x))^2),y-b*sin(2*pi*x);");
		f[i]->setParameter("pi",interval::pi());
		fd[i]=new DMap("par:pi,a,b;var:x,y;fun:x+a*(1-(y-b*sin(2*pi*x))^2),y-b*sin(2*pi*x);");
		fd[i]->setParameter("pi",4.0*atan(1.0));
		counter[i]=0;
	}

	int i;
	int done=0;
	cout << "progress (the program will finish at 1.0): " << endl;
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
				counter[id]++;
			}
		}
		if(Na/100 > 0)
		{
			if( (i % Na/100) == 0)
			{
				done++;
				cout << done/100.0 << endl;
			}
		}else{
			done++;
			cout << done/double(Na) << endl;
		}
	}
	int total=0;
	for(int i=0;i<N_of_threads;i++) total=total+counter[i];
	cout << "Total number of parameter pairs for which we have unbounded diffusion in NTSF: " << total << endl << endl; 
}

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

void ProofOfTheorem_1_2()
{
	auto start = chrono::high_resolution_clock::now(); // Start time
	
	validateChaosInConservativeMaps();

	auto end = chrono::high_resolution_clock::now(); // End time
	chrono::duration<double> duration = end - start;
  	cout << "Execution time: " << duration.count() << " seconds" << endl;
}

void ProofOfTheorem_1_3(int N)
{
	auto start = chrono::high_resolution_clock::now(); // Start time
	
	validateDiffusionNoTwistSF(N);

	auto end = chrono::high_resolution_clock::now(); // End time
  	chrono::duration<double> duration = end - start;
  	cout << "Execution time: " << duration.count() << " seconds" << endl;
}

void ProofOfTheorem_1_4(int accuracy)
{
	auto start = chrono::high_resolution_clock::now(); // Start time
	
	validateChaosInStandardMapCluster(accuracy);

	auto end = chrono::high_resolution_clock::now(); // End time
  	chrono::duration<double> duration = end - start;
  	cout << "Execution time: " << duration.count() << " seconds" << endl;
}

int main(int argc, char* argv[])
{	
	cout.precision(10);
	try
	{	
		if((argc==1) or (argc==2 and std::atoi(argv[1])==0))
		{ 
			ProofOfTheorem_1_2();
		}
		
		if(argc==3) 
  		{
  			int whichProof = std::atoi(argv[1]);
  			int n = std::atoi(argv[2]);

  			if(whichProof==1)
  			{
  				if(n<0) return 1;
  				ProofOfTheorem_1_3(n);
  			}
  			if(whichProof==2)
  			{
  				if(n<0) return 1;
  				if(n>3) return 1;
  				ProofOfTheorem_1_4(n);
  			}
  		}
		
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  	return 0;
} 
