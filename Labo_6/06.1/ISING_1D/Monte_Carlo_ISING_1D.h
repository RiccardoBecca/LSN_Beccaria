

#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
#include <cstdio>
#include <vector>
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig, iu2;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g, stima_u2;
double err_u,err_c,err_m,err_x,err_g;

//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk, metro, restart;

//functions
void Input(int);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);
double Min(double,double);
double Error(std::vector<double>, std::vector<double>, int);

#endif


