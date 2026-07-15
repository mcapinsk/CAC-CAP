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

	// Removing (author's) results which are downloaded from github:
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

int max(vector<int> a)
{
	int k=a[0];
	for(int i=1;i<(int)a.size();i++)
	{
		if(k<a[i]) k=a[i];
	}
	return k;
}

interval validateChaosInStandardMapCluster(int Debth)
{
	int nA=700;
	int nB=100;

	interval A(3.0,10.0);
	interval B(0.01,0.8);
	ofstream plotFile("thm_1_4/dissipative-area.txt");

	int N_of_threads=omp_get_max_threads();
	cout << "Number of threads: " << N_of_threads << endl;

	vector<IStandardMap*> F(N_of_threads);
	vector<DStandardMap*> Fd(N_of_threads);
	vector<chaosProofParameters> par(N_of_threads);
	vector<interval> Area(N_of_threads),a(N_of_threads),dArea(N_of_threads);
	vector<int> NofIterates(N_of_threads);

	for(int i=0;i<N_of_threads;i++)
	{
		F[i]=new IStandardMap();
		Fd[i]=new DStandardMap();
		par[i].B = interval(0); // level to go above/below (for below we take -B)
		                        // The B will be updated during the proof, since it depends on the parameters of the map.
		par[i].maxIteratesUpDown=40;

		// This will count the maximum length of trajectory 
		// needed to go above B or below -B.
		NofIterates[i]=0;
	}
	
	cout << "progress (the program will finish at 1.0): " << endl;

	vector<interval> BArea(nB);
	
	for(int j=0;j<nB;j++)
	{
		interval b=part(B,nB,j);
		BArea[j]=0;
		
		vector<interval> bAreaPart(N_of_threads);

		int i;
		#pragma omp parallel for private(i)
		for(i=0;i<nA;i++)
		{
			int id=omp_get_thread_num();
			a[id]=part(A,nA,i);
			dArea[id] = chaoticArea(*(F[id]),*(Fd[id]),a[id],b,par[id],Debth,NofIterates[id]);
			bAreaPart[id] = bAreaPart[id] + dArea[id];
			Area[id] = Area[id] + dArea[id];
		}

		for(int k=0;k<N_of_threads;k++) BArea[j]=BArea[j]+bAreaPart[k];
		plotFile << b.leftBound() << " " <<  (BArea[j]/((b.right()-b.left())*(A.right()-A.left()))).leftBound() << endl;
		plotFile << b.rightBound() << " " << (BArea[j]/((b.right()-b.left())*(A.right()-A.left()))).leftBound() << endl;
		cout << j << " out of " << nB << endl;
	}

	interval total_Area(0);
	for(int k=0;k<N_of_threads;k++) total_Area = total_Area+ Area[k];
	interval total_Area2(0);
	for(int k=0;k<nB;k++) total_Area2 = total_Area2+ BArea[k];
	cout << "total area validated        : " << total_Area << endl;
	cout << "total area validated        : " << total_Area2 << endl;
	cout << "percentage of area validated: " << total_Area/((A.right()-A.left())*(B.right()-B.left())) << endl;
	cout << "largest number of iterates used: " << max(NofIterates) << endl;

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
  				if(n>2) return 1;
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
