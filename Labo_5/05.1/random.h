

#ifndef __Random__
#define __Random__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "funzionebase.h"

class Random {

private:
    int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
    // constructors
    Random();
    // destructor
    ~Random();
    // methods
    void SetRandom(int * , int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
    double Exp(double lamda);
    double Lorentz(double mu, double gamma);
    double Pi();
    double d_prob(double a, double b, double f_max, FunzioneBase* f);
    
    void Passo(vector<double> & v);
    void Passo_continuo(vector<double> & v);
};



#endif // __Random__

