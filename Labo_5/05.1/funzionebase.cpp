#include "funzionebase.h"
#include <math.h>

using namespace std;

Coseno::Coseno(){
    m_a=1;
    m_b=1;
    m_c=0;
}

Coseno::Coseno(double a, double b, double c){
    m_a=a;
    m_b=b;
    m_c=c;
}

Parabola::Parabola(){
    m_a=1;
    m_b=1;
    m_c=1;
}

Parabola::Parabola(double a, double b, double c){
    m_a=a;
    m_b=b;
    m_c=c;
}

Retta::Retta(){
    m_m=1.;
    m_q=0.;
}

Retta::Retta(double m, double q){
    m_m=m;
    m_q=q;
}

H_gs_density::H_gs_density(){
    m_bohr_radius= 0.0529*pow(10,-9);
}

H_gs_density_2::H_gs_density_2(){
    m_bohr_radius= 0.0529*pow(10,-9);
}

double H_gs_density_2::Eval(vector<double> v) const{
    double theta;
    if(v[0]==0 && v[1]>0){
        theta= M_PI*0.5;
    }
    if(v[0]==0 && v[1]<0){
        theta= 3*M_PI*0.5;
    }
    if(v[0]>0 && v[1]>=0){
        theta=atan(v[1]/v[0]);
    }
    if(v[0]>0 && v[1]<0){
        theta=atan(v[1]/v[0]) + 2*M_PI;
    }
    if(v[0]<0 && v[1]>0){
        theta=atan(v[1]/v[0]) + 2*M_PI;
    }
    if(v[0]<0 && v[1]<=0){
        theta=atan(v[1]/v[0]) +M_PI;
    }
    
    return pow(m_bohr_radius, -5)*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*exp(-(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/m_bohr_radius))*cos(theta)*cos(theta)/(32*M_PI);
    
}
