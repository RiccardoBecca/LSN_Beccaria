#ifndef __NVT__
#define __NVT__

//#include "random.h"
#include <vector>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> // for sort()
#include <math.h>
#include <cstdlib>
#include <random>
#include <cmath>

using namespace std;


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
};


int PBC(int,int);
double Min(double,double);

class Cities{
    vector<vector<double>> m_city;
    int m_numero;//32 nel nostro caso
public:
    Cities();
    Cities(int,int,Random); 
    void Print(void);
    Cities(vector<vector<double>> v);
    double Get(int a, int b){return m_city[a][b];}
    int Getlen(){return m_numero;}
    
};

class Figlio{
    vector<int> m_son;
    double m_len;
    
public:
    Figlio(int);
    int Get(int a){return m_son[a];}
    double Getlen(){return m_len;}
    void Print(void);
    void Check(int);
    void Copy(Figlio);
    int Size(void){return m_son.size();}
    void Measure_len(Cities);
    void Set(int a, int b){m_son[a]=b;}
    void Setlen(double a){m_len=a;}
    void Swap(int a, int b){swap(m_son[a],m_son[b]);}
};

class Generazione{
    vector<Figlio> m_gen;
    int m_lunghezza; //32 nel nostro caso
    int m_num_figlio; //100 nel mio caso
    double m_prob_cross, m_prob_m1, m_prob_m2, m_prob_m3, m_prob_m4;
    Random m_rand;
    int m_T;
    int m_accepted;
public:
    Generazione(int,int,double,double,double,double,double,Random,int);
    double Get_Best_Len();
    double Get_Mean_Len();
    void Print(void);
    void Check(void);
    void Measure(Cities);
    void Sorting();
    void Crossing(Cities, Generazione);
    void Mutation(Cities);
    void Mutation_swap(Cities);
    void Mutation_permut(Cities);
    void Mutation_shift(Cities);
    void Mutation_inversion(Cities);
    Figlio Get_son(int a){return m_gen[a];}
    void Copy(Generazione);
    void Metro(Cities);
    void Metro_swap(Cities);
    void Metro_permut(Cities);
    void Metro_shift(Cities);
    void Metro_inversion(Cities);
    int Get_accepted(void){return m_accepted;}
    void Set_accepted(int a){m_accepted=a;}
    void Set_T(double T){m_T=T;}
};
bool sortcol(Figlio & v1, Figlio & v2);




#endif

