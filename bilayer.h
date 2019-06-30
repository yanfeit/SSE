#ifndef __bilayer__
#define __bilayer__

#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <string>
#include <sstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;

typedef vector<vector<int> > Matrix;
typedef vector<vector<double> > dMatrix;

struct Model{
	int lx; int ly;
	double beta; int mpps;
	int warmsteps; int mcsteps;
	double t_parallel; double t_perp;
	double interaction; double chempot;
	int seed;
	
	// derived member from input file
	int ns; double c; // const c, keep in mind to add N_b*c to total energy at last
	int nb;
	
	Matrix bsites; vector<double> weight;
	vector<double> table;
};

struct Configuration{
	
	int msl; int nh;
	vector<int> state;
	vector<int> opstring;
	Matrix leg;
	vector<int> vertexlist;
	vector<int> frststateop;
	vector<int> laststateop;
	
};

struct Observables{
	vector<double> density;
	vector<double> energy;
	vector<double> rhosTT;
	vector<double> rhosTB;
	vector<double> rhos;
		
};

void readfile(Model &, int, char **);

void makelattice(Model &);

void calweight(Model &);

dMatrix solution(const vector<double> &w);

void askspace(const Model &, Configuration &);

void diagupdate(const Model &, Configuration &);

void adjustcutoff(const Model &, Configuration &);

void linkvertices(const Model &, Configuration &);

void caltable(Model &);

// some function used to compute the probablity table
// ---------------------------------------------------------------------------------//-
double bounce(const Model &, int kind, int type, int np, int nq, int crement);      //-
double reverse(const Model &, int kind, int type, int np, int nq, int crement);     //-
double straight(const Model &, int kind, int type, int np, int nq, int crement);	//-
double jump(const Model &, int kind, int type, int np, int nq, int crement);	    //-
//----------------------------------------------------------------------------------//-

int determinetype(const int, const int);

int rotation(int x, int y);

int loopupdate(const Model &, Configuration &);

void updatestate(const Model &, Configuration &);

void initialization(Model &, Configuration &, int, char **);

void warmup(const Model &, Configuration &);

int computelooptimes(const Model &, Configuration &);

void display(const Model &, const Configuration &);

void simulation(const Model &, Configuration &, Observables &, int);

double computedensity(const Configuration &);

double computeenergy(const Model &, const Configuration &);

vector<double> computewind(const Model &, const Configuration &);

void writetofile(const Model &, const Observables &);

// functions used to do data analysis

vector<double> jackknife(const vector<double> &);

vector<double> binning(const vector<double> &);

vector<double> computek(const Model &, const vector<double> &);

#endif