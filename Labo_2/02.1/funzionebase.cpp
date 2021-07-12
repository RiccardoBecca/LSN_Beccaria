#include "funzionebase.h"


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
