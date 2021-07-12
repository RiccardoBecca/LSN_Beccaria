#ifndef __lib__
#define __lib__

#include <vector>
#include "funzionebase.h"
#include "random.h"


double Error(std::vector<double>, std::vector<double>, int n);

double Error(double a, double b, int n);

double max(double a, double b);
double min(double a, double b);

void Metropolis_unif(std::vector<double> &v, FunzioneBase* f, Random* rnd, double s, int & n_hit, int &n_tot);

void Metropolis_gauss(std::vector<double> &v, FunzioneBase* f, Random* rnd, double s, int & n_hit, int &n_tot);

#endif
