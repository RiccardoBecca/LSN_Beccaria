
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "lib.h"
#include "funzionebase.h"
#include "integral.h"

using namespace std;
 
int main (int argc, char *argv[]){


//    Random*rnd;
    Random* rnd=new Random();
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
//            rnd.SetRandom(seed,p1,p2);
             rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


    
   
    ofstream myfile;
    Coseno* mycos = new Coseno(M_PI*0.5, M_PI*0.5, 0.);
    Integral* myintegral= new Integral(0.,1., mycos, rnd);

    int M=pow(10,4);
    int N=300;
    int L=M/N;
    double value=0.;
    double value2=0.;
    double error=0.;
    double x=0;
    
    
    myfile.open("integral.txt");
    for(int i=0; i<N; i++){
        x=myintegral->Ave(L);
        value=(value*i+x)/double(i+1);
        value2=(value2*i+x*x)/double(i+1);
        error= Error(value, value2, i);
        myfile << value << " " << value2 << " " << error <<endl;
    }

    myfile.close();

    myfile.open("imp_samp_retta.txt");
    Retta* myretta = new Retta(-2., 2.);
    for(int i=0; i<N; i++){
        x=myintegral->Imp_samp(L, 2. , myretta);
        value=(value*i+x)/double(i+1);
        value2=(value2*i+x*x)/double(i+1);
        error= Error(value, value2, i);
        myfile << value << " " << value2 << " " << error <<endl;
    }
    myfile.close();


    myfile.open("imp_samp_parab.txt");
    Parabola* myparab = new Parabola(-1.,0.,4./3.);
    for(int i=0; i<N; i++){
        x=myintegral->Imp_samp(L, 4./3. ,myparab);
        value=(value*i+x)/double(i+1);
        value2=(value2*i+x*x)/double(i+1);
        error= Error(value, value2, i);
        myfile << value << " " << value2 << " " << error <<endl;
    }
    myfile.close();

    
    
    
    rnd->SaveSeed();
    return 0;
}


