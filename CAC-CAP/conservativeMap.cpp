#include "conservativeMap.h"
#include "validatedTrajectory.h"
#include "utils.h"

// This routine returns the formuale of conservative maps considered
// in the paper.
string mapFormula(int i)
{
	// Below formulae were generated in:
	//   Formulae_conservative_map/MathematicaFormulaeComputation.nb

	// ==========Table 1==========
	if(i==0) return "x+y,y+sin(2*pi*(x+y))+c";
	//========
	if(i==1) return "x+y,y+(1-sin(2*pi*(x+y)))*sin(2*pi*(x+y))+c";
	//========
	if(i==2) return "x+y,y+sin(sin(2*pi*(x+y)))/cos(sin(2*pi*(x+y)))+c";
	//========
	if(i==3) return "x+y,y+3*log(2+sin(2*pi*(x+y)))+c";
	//========
	if(i==4) return "x+y,-1+y+exp(sin(2*pi*(x+y)))+c";
	
	//==========Table 2==========
	if(i==5) return "x+sin(2*pi*y),y+sin(2*pi*(x+sin(2*pi*y)))+c";
	//========
	if(i==6) return "x+sin(2*pi*y),y+(1-sin(2*pi*(x+sin(2*pi*y))))*sin(2*pi*(x+sin(2*pi*y)))+c";
	//========
	if(i==7) return "x+sin(2*pi*y),y+sin(sin(2*pi*(x+sin(2*pi*y))))/cos(sin(2*pi*(x+sin(2*pi*y))))+c";
	//========
	if(i==8) return "x+sin(2*pi*y),y+3*log(2+sin(2*pi*(x+sin(2*pi*y))))+c";
	//========
	// i==9
	return "x+sin(2*pi*y),-1+y+exp(sin(2*pi*(x+sin(2*pi*y))))+c";
}

// This routine computes the integral of v(sin(2*pi*x))−c.
// It is used to find the constant c for the considered function v.
interval integral(int i)
{
	// in these cases c=0
	if((i==0) or (i==2) or (i==5) or (i==7)) return interval(0);

	// in these cases c=1/2
	if((i==1) or (i==6)) return -interval(1)/interval(2);

	IMap w;
	if((i==3) or (i==8)) w=IMap("par:pi;var:x;fun:3*log(2+sin(2*pi*x));");
	if((i==4) or (i==9)) w=IMap("par:pi;var:x;fun:exp(sin(2*pi*(x)))-1;");
	w.setParameter("pi",interval::pi());
	int N=100000;
	interval Integral(0);
	for(int i=0;i<N;i++)
	{
		interval x(0,1);
		x=part(x,N,i);
		Integral=Integral+w(IVector({x}))[0]*(x.right()-x.left());
	}
	return Integral;
}

/////////////////////////////////
// This function considers N points on y=-B and finds a 
// shortest trajectory going above B.
//
// Note: the function always returns some sequence of points,
// but if the length is equal to max_n, then such
// trajectory does not go above B. This is not a problem,
// since each trajectory is additionally validated in 
// validateChaosInConservativeMaps().
vector<DVector> findShortestTrajectoryUp(DMap &f,double B,int N,int max_n)
{
	int shortest_path=max_n;
	int best_i=0;
	for(int i=0;i<N;i++)
	{
		DVector x({(double)i/(double)N,-B});
		int j=0;
		while(j<shortest_path)
		{
			x=f(x);
			if(x[1]>B+0.01)
			{
				shortest_path=j+1;
				best_i=i;
			}
			j++;
		}
	}
	vector<DVector> q;
	DVector x({(double)best_i/(double)N,-B});
	q.push_back(x);
	for(int j=0;j<shortest_path;j++)
	{
		x=f(x);
		q.push_back(x);
	}
	return q;
}

// We choose default N=10000 and n_max=100 to search for 
// trajectories which go up in validateChaosInConservativeMaps().
vector<DVector> findShortestTrajectoryUp(DMap &f,double B)
{
	return findShortestTrajectoryUp(f,B,10000,100);
}

// This is a technical function, used just for plots:
void plotTrajectory(vector<DVector> q,int i)
{
	ofstream file("plots/conservative"+to_string(i)+".txt");
	for(int k=0;k<q.size();k++) file << q[k][0] << " " << q[k][1] << endl;
	file.close();
}

// This is the theorem in which we validate Theorem 1.2.
void validateChaosInConservativeMaps()
{
	for(int i=0;i<10;i++)
	{
		if(i==0) cout << "======== Table 1 ========" << endl;
		if(i==5) cout << "======== Table 2 ========" << endl;
		interval Int=integral(i); // this is -c

		// The DMap (double map) is used to find a candidate
		// for the trajectory.
		DMap f("par:pi,c;var:x,y;fun:"+mapFormula(i)+";");
		f.setParameter("pi",4.0*atan(1.0));
		f.setParameter("c",-Int.mid().leftBound());

		// We take B=5
		interval B(5.0);

		// This is the candidate for a trajectory going from
		// y=-B to above y=B.
		vector<DVector> q=findShortestTrajectoryUp(f,B.leftBound());

		plotTrajectory(q,i);
		bool result;
		
		vector<IVector> q_validated;
		if(q.size()<100)
		{
			IVector q0(2);
			q0[0]=q[0][0];
			q0[1]=-B;
			IVector v0({interval(1),interval(0)});

			IMap F("par:pi,c;var:x,y;fun:"+mapFormula(i)+";");
			F.setParameter("pi",interval::pi());
			F.setParameter("c",-Int);
			interval h0=interval(0.0);

			result=validateTrajectoryUp(q,q_validated,B,F,q0,v0,h0);
			// if the trajectory is validated the result is set to 1.
		}else
		{
			// the proof failed.
			result=0;
		}
		cout << i % 5 +1 << ". Chaos in the conservative map: (" << mapFormula(i) << ") validated." << endl;
		if(result==1)
		{
			cout << "c                   : " << -Int << endl;
			cout << "Initial point       : {" << q_validated[0][0].mid().leftBound() << "+" << q_validated[0][0]-q_validated[0][0].mid() <<","<<q_validated[0][1] << "}"  << endl;
			cout << "Length of trajectory: " << q.size()-1 << endl << endl; 
		}else
		{
			cout << "NOT validated." << endl;
		}
	}
}
