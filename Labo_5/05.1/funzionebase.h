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


class H_gs_density: public FunzioneBase {
public:
    H_gs_density();
    ~H_gs_density() {;};
    
    virtual double Eval(vector<double> v) const {return pow(m_bohr_radius, -3)*exp(-2*(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))/m_bohr_radius)/M_PI;}
    virtual double Eval (double x) const {
        cout << "This function needs a vector, not a double. There is something wrong" << endl;
        return 0;
    };

    void Set_bohr_radius (double a){m_bohr_radius=a;}
    double Get_bohr_radius (){return m_bohr_radius;}
    
private:
    double m_bohr_radius;
};

class H_gs_density_2: public FunzioneBase {
public:
    H_gs_density_2();
    ~H_gs_density_2() {;};
    
    //virtual double Eval(vector<double> v) const {return pow(m_bohr_radius, -3)*exp(-2*(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))/m_bohr_radius)/M_PI;}
    virtual double Eval(vector<double> v) const;
    virtual double Eval (double x) const {
        cout << "This function needs a vector, not a double. There is something wrong" << endl;
        return 0;
    };

    void Set_bohr_radius (double a){m_bohr_radius=a;}
    double Get_bohr_radius (){return m_bohr_radius;}
    
private:
    double m_bohr_radius;
};

#endif
