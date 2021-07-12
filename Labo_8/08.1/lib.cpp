
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "lib.h"
#include "random.h"
#include "funzionebase.h"

using namespace std;


double Error(std::vector<double> ave, std::vector<double> av2, int n){
    if(n==0){
        return 0;
    }
    else{
        return sqrt((av2[n]-ave[n]*ave[n])/double(n));
    }
}

double Error(double ave, double av2, int n){
    return sqrt((av2-ave*ave)/double(n));
}

double max(double a, double b){
    if (a<b){
        return b;
    }
    else{
        return a;
    }
}

double min(double a, double b){
    if (a<b){
        return a;
    }
    else{
        return b;
    }
}

void Metropolis_unif(vector<double> &v, FunzioneBase* f, Random* rnd, double s, int & n_hit, int &n_tot){
    vector<double> w (3,0);
    w[0]= v[0]+s*rnd->Rannyu(-1,1);
    w[1]= v[1]+s*rnd->Rannyu(-1,1);
    w[2]= v[2]+s*rnd->Rannyu(-1,1);
    
    double alpha=min(1, (f->Eval(w))/(f->Eval(v)));
    
    double r=rnd->Rannyu(0,1);
    if(r<alpha){
        v[0]=w[0];
        v[1]=w[1];
        v[2]=w[2];
        n_hit+=1;
    }
    else{
        v[0]=v[0];
        v[1]=v[1];
        v[2]=v[2];
    }
    n_tot +=1;
}

void Metropolis_gauss(vector<double> &v, FunzioneBase* f, Random* rnd, double s, int & n_hit, int &n_tot){
    vector<double> w (3,0);
    w[0]= rnd->Gauss(v[0],s);
    w[1]= rnd->Gauss(v[1],s);
    w[2]= rnd->Gauss(v[2],s);
    
    double alpha=min(1, (f->Eval(w))/(f->Eval(v)));
    
    double r=rnd->Rannyu(0,1);
    if(r<alpha){
        v[0]=w[0];
        v[1]=w[1];
        v[2]=w[2];
        n_hit+=1;
    }
    else{
        v[0]=v[0];
        v[1]=v[1];
        v[2]=v[2];
    }
    n_tot +=1;
}


double Metropolis_psi(double x, FunzioneBase* f, Random* rnd, int & n_hit, int &n_tot){
    
    double x2;
    x2 = x+rnd->Rannyu(-1,1);
    n_tot +=1;
    double alpha=min(1, (f->Eval(x2))*(f->Eval(x2))/((f->Eval(x))*(f->Eval(x))));
    double r=rnd->Rannyu(0,1);
    
    if(r<= alpha){
        n_hit+=1;
        return x2;
    }
    else{
        return x;
    }
}
