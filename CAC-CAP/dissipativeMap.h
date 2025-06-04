#ifndef dissipativeMap_h
#define dissipativeMap_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

// We create a class to store a dissipative map.
// This class essentially mimics the functionality of DMap from CAPD,
// that is, it is used to compute the image and derivative of a function
//    R^n -> R^n.
// The main reason for creating such class is that for our problem a map
// is associated with two fixed points, which are used in the proof.
// We have therefore added a virtual function 
//    DVector FixedPoint(int i)
// (We have two fixed points for i=0,1)
// This is an abstract class, which can be used create subclasses for
// particular types of maps. We have created just one subclass, the dissipative SM,
// but the method can be applied to other types of maps by creating a subclass 
// and performing the proof with it.
class DDissipativeMap
{
public:
	DMap *map;

	DDissipativeMap(string formula){map = new DMap(formula);}

	void setParameter(string par,double var){map->setParameter(par,var);}

	DVector value(const DVector &x) const {return (*map)(x);}
	DMatrix derivative(const DVector &x) const {return (*map)[x];}

	DVector operator()(const DVector &x) const {return (*map)(x);}
	DMatrix operator[](const DVector &x) const {return (*map)[x];}

	virtual DVector FixedPoint(int i)=0;
	virtual ~DDissipativeMap(){delete map;}
};

// In our dissipative SM we have a family of pairs of fixed points which 
// depend on the parameter kappa. We therefore add a parameter to the class.
// We can change that parameter and use a different pair of fixed point in
// the proof.
class DStandardMap: public DDissipativeMap
{
private:
	double PI;
	int kappa; // choice of fixed points pair
	double a,b;
public:
	DStandardMap(double a,double b);
	DStandardMap();

	void set_kappa(int k){kappa=k;}
	void set_a(double a);
	void set_b(double b);

	int get_kappa(){return kappa;}
	double get_a(){return a;}
	double get_b(){return b;}

	DVector FixedPoint(int i);
};

////////////////////////////////////////////////////

// This class is an interval arithmetic version of the class DDissipativeMap.
// It has the same functionality.
class IDissipativeMap
{
public:
	IMap *map;

	IDissipativeMap(string formula){map = new IMap(formula);}

	void setParameter(string par,interval var){map->setParameter(par,var);}

	IVector value(const IVector &x) const {return (*map)(x);}
	IMatrix derivative(const IVector &x) const {return (*map)[x];}

	IVector operator()(const IVector &x) const {return (*map)(x);}
	IMatrix operator[](const IVector &x) const {return (*map)[x];}

	virtual ~IDissipativeMap(){delete map;}

	virtual IVector FixedPoint(int i)=0;
};

// This class is an interval arithmetic version of the class DStandardMap.
// It has the same functionality.
class IStandardMap: public IDissipativeMap
{
private:
	interval PI;
	int kappa; // choice of fixed points pair
	interval a,b;
public:
	IStandardMap(interval a,interval b);
	IStandardMap();

	void set_kappa(int k){kappa=k;}
	void set_a(interval a);
	void set_b(interval b);

	int get_kappa(){return kappa;}
	interval get_a(){return a;}
	interval get_b(){return b;}

	IVector FixedPoint(int i);
};

/////////////////////////////////////////

// This map is used to compute the derivative in clocal coordinates.
// Its main purpose is for the validation of cone conditions.
// We consider a map f:R^n -> R^n, an n x n matrix A, and two points x0,x1. 
// The map f in local coordinates is given as
//     x-> Ainv*(f(x0+A*x)-x1).
// The class is used to compute the image and the derivative of this map in 
// the local coordinates.
class localMap
{
private:
	IMatrix A,Ainv;
	IVector x0,x1;
	IDissipativeMap *f;
public:
	IMatrix get_A(){return A;}
	IVector image(IVector x){return Ainv*((*f)(x0+A*x)-x1);}
	IMatrix derivative(IVector x){return Ainv*(*f)[x0+A*x]*A;}

	localMap(IVector x0,IVector x1,IMatrix A,IMatrix Ainv,IDissipativeMap &f);
	
	IVector operator()(IVector x){return image(x);}
	IMatrix operator[](IVector x){return derivative(x);}
};

/////////////////////////////////////////

// Ihe parallel shooting Krawczyk operator depend on a number of parameters.
// We bundle them in a single class. This class serves as a container of these 
// parameters.
class chaosProofParameters
{
public:
	interval B; // level to go above/below (for below we take -B)
	interval r; // how far from the fixed point we use the cone
	interval L; // initial cone slope
	double rho; // where position the curve at the end of the cone
	int maxIteratesUpDown; // the maximal length of trajectory we search for
};

//////////////////////////////////////////
// This function is of main interest. It takes a small box of parameters, which
// is a cartesian product of the two intervals a and b. The function tries to validate
// the appropriate conditions, namely the existence of trajectories which go
// up/down fromt the fixed points to above/below the level L, and to validate 
// the segment conditions. If the map fails, then it subdivides the intervals a,b
// and tries again. If for some parameters it succeeds, it counts them towards
// the validated area. 
// - The number of recursive subdivisions is chosen by the variable Debth.
// - The map F is used to validate the needed trajectories and segments
// - The map Fd is used to perform non-rigorous simulations to find the
//   candidates for the appropriate trajectories to be validated.
interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Debth);

#endif