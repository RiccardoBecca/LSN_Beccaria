#ifndef __Integral__
#define __Integral__

#include <cmath>
#include <iostream>
#include <cstdlib>
#include "funzionebase.h"
#include "random.h"


using namespace std;

class Integral{
public:
    Integral(double a, double b, FunzioneBase* f, Random* rnd);
    
    double Ave(unsigned int punti);
    double Imp_samp(unsigned int punti, double d_max, FunzioneBase* d); // voglio solo dargli in input la densità di prob con cui volgio campionare
    
private:
    double m_min;
    double m_max;
    double m_integral;
    int m_sign;
    FunzioneBase * m_f;
    Random* m_rnd;
    

};
    
#endif
