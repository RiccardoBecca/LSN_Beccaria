#ifndef __FunzioneBase__
#define __FunzioneBase__

#include <cmath>
#include <iostream>
#include <cstdlib>


using namespace std;

class FunzioneBase {

 public:
    virtual double Eval (double x) const=0;
};

class Coseno: public FunzioneBase {
public:
    Coseno();
    Coseno(double a, double b, double c);
    ~Coseno() {;};
    
    virtual double Eval(double x) const {return m_a*cos(m_b*x+m_c);}

    void SetA (double a){m_a=a;}
    void SetB (double b){m_b=b;}
    void SetC (double c){m_c=c;}
    double GetA (){return m_a;}
    double GetB (){return m_b;}
    double GetC (){return m_c;}
    
private:
    double m_a;
    double m_b;
    double m_c;
};

class Parabola: public FunzioneBase {
public:
    Parabola();
    Parabola(double a, double b, double c);
    ~Parabola() {;};
    
    virtual double Eval(double x) const {return m_a*x*x+m_b*x+m_c;}
    
private:
    double m_a;
    double m_b;
    double m_c;
    
    
};

class Retta: public FunzioneBase {
public:
    Retta();
    Retta(double m, double q);
    ~Retta() {;};
    
    virtual double Eval(double x) const {return m_m*x+m_q;}
    
private:
    double m_m;
    double m_q;
};



#endif
