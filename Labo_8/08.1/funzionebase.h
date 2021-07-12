#ifndef __FunzioneBase__
#define __FunzioneBase__

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>


using namespace std;

class FunzioneBase {

 public:
    virtual double Eval (double x) const=0;
    virtual double Eval (vector<double> v) const=0;
    
    virtual double Second_derivate(double x) const=0;
};

class Coseno: public FunzioneBase {
public:
    Coseno();
    Coseno(double a, double b, double c);
    ~Coseno() {;};
    
    virtual double Eval(double x) const {return m_a*cos(m_b*x+m_c);}
    virtual double Second_derivate(double x) const {return 0;}

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
    virtual double Second_derivate(double x) const {return 0;}
    
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
    virtual double Second_derivate(double x) const {return 0;}
    
private:
    double m_m;
    double m_q;
};


class Psi_T: public FunzioneBase {
public:
    Psi_T();
    Psi_T(double s, double m);
    ~Psi_T() {;};
    
    virtual double Eval(double x) const {return exp(-pow((x-m_mu),2)/(2*m_sigma*m_sigma) )+ exp(-pow((x+m_mu),2)/(2*m_sigma*m_sigma));}
    virtual double Eval(vector<double> v) const {return 0;}
    virtual double Second_derivate(double x) const {return -(exp(-pow(x-m_mu,2)/(2*m_sigma*m_sigma)) + exp(-pow(x+m_mu,2)/(2*m_sigma*m_sigma)) )/(m_sigma*m_sigma) +( pow(x-m_mu,2)*exp(-pow(x-m_mu,2)/(2*m_sigma*m_sigma)) +pow(x+m_mu,2)*exp(-pow(x+m_mu,2)/(2*m_sigma*m_sigma)) )/(pow(m_sigma,4));}


    void Set_m_sigma (double s){m_sigma=s;}
    double Get_m_sigma (){return m_sigma;}
    void Set_m_mu (double m){m_sigma=m;}
    double Get_m_mu (){return m_mu;}
    
private:
    double m_sigma;
    double m_mu;
};


class Double_well: public FunzioneBase {
public:
    Double_well();
    Double_well(double a, double b);
    ~Double_well() {;};
    
    virtual double Eval(double x) const {return (m_b*pow(x,2)+m_a*pow(x,4)) ;}
    virtual double Eval(vector<double> v) const {return 0;}
    virtual double Second_derivate(double x) const {return 12*m_a*x*x+2*m_b;}

    void SetA (double a){m_a=a;}
    void SetB (double b){m_b=b;}
    double GetA (){return m_a;}
    double GetB (){return m_b;}
    
private:
    double m_a;
    double m_b;
};

#endif
