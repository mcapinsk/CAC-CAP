#ifndef dissipativeMap_h
#define dissipativeMap_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

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

class IStandardMap: public IDissipativeMap
{
private:
	interval PI;
	int kappa; // choice of fixed points pair
	interval a,b;
public:
	IStandardMap(interval a,interval b);
	IStandardMap();
	//IStandardMap& operator=(const IStandardMap& other);

	void set_kappa(int k){kappa=k;}
	void set_a(interval a);
	void set_b(interval b);

	int get_kappa(){return kappa;}
	interval get_a(){return a;}
	interval get_b(){return b;}

	IVector FixedPoint(int i);
};

/////////////////////////////////////////

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

class chaosProofParameters
{
public:
	interval B; // level to go above/below (for below we take -B)
	interval r; // how far from the fixed point we use the cone
	interval L; // initial cone slope
	double rho; // where position the curve at the end of the cone
	int maxIteratesUpDown;
};

//////////////////////////////////////////

bool chaosInStandardMapNonRigorous(DStandardMap &F,int maxIteratesUpDown,double r,int &NofIterates);

interval chaoticArea(IStandardMap &F,DStandardMap &Fd,interval a,interval b,chaosProofParameters &par,int Debth);

#endif