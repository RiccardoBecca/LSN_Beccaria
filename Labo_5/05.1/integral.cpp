#include "integral.h"


#include <algorithm>

using namespace std;

Integral::Integral(double a, double b, FunzioneBase* f, Random* rnd){
    m_min=min(a,b);
    m_max=max(a,b);
    m_f=f;
    m_rnd=rnd;
    if(a>b) m_sign=-1.;
    else m_sign=1.;
    m_integral=0.;
}

double Integral::Ave(unsigned int punti){
    m_integral = 0.;
    for (unsigned int i=0; i<punti; i++) {
        m_integral += m_f->Eval(m_rnd->Rannyu(m_min, m_max));
    }
    m_integral = (double)m_sign*(m_max-m_min)*m_integral/(double)punti;
    return m_integral;
}


double Integral::Imp_samp(unsigned int punti, double d_max, FunzioneBase* d){
    m_integral = 0.;
    for (unsigned int i=0; i<punti; i++) {
        double x=m_rnd->d_prob(m_min,m_max, d_max, d);
        m_integral += (m_f->Eval(x))/(d->Eval(x));
    }
    m_integral = (double)m_sign*(m_max-m_min)*m_integral/(double)punti;
    return m_integral;
}
