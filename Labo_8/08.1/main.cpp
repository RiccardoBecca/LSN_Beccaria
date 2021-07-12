
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include "random.h"
#include "lib.h"
#include "funzionebase.h"
#include "integral.h"

using namespace std;
 
int main (int argc, char *argv[]){


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
             rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    
    ofstream myfile;
    myfile.open("Best_Psi.txt");
    double x=0;
    Double_well* myPot = new Double_well(1., -2.5);
    int nhit = 0;
    int ntot = 0;
    

    int L = 10000;
    int N = 100;
    double integral=0;
    double minus=100;
    double s_min=0;
    double m_min=0;


    for(double m=0.78; m<0.83; m+=0.001){
        for(double s=0.6; s<0.65; s+=0.001){
            
            double value =0;
            double value2=0;
            double error =0;
            nhit = 0;
            ntot = 0;
            
            Psi_T* myPsi  = new Psi_T(s, m);
            for(int k=0; k<1000; k++){ //equilibration
                x = Metropolis_psi(x, myPsi, rnd, nhit, ntot);
            }
            for(int n=0; n<N; n++){
                integral=0;
                for(int i=0; i<L; i++){
                    x = Metropolis_psi(x, myPsi, rnd, nhit, ntot);
                    integral += (-0.5*(myPsi->Second_derivate(x)) + (myPot->Eval(x))*(myPsi->Eval(x)) ) /(myPsi->Eval(x));
                }
                integral=integral/double(L);
                value=(value*n+integral)/double(n+1);
                value2=(value2*n+integral*integral)/double(n+1);
                error= Error(value, value2, n);
                if(n==99){
                    myfile << s << " " << m << " " << value << endl;
                    cout << "Sigma = " << s << "   Mu = " << m << "     Integral = " << value << endl;
                    if(value<minus){
                        minus=value;
                        s_min=s;
                        m_min=m;
                    }
                }
                
            }
        }
    }
    cout << "=====================" << endl << "GS EIGENVALUE = " << minus << endl <<"Sigma = " << s_min << endl << "Mu = " << m_min << endl << "=====================" << endl;
    
    myfile.close();
    

    
    
    double xmax=2.75;
    
    ofstream gs_eigenvalue;
    gs_eigenvalue.open("Ground_state_eigenvalue.txt");
    myfile.open("Ground_state_eigenstate.txt");
    
    Psi_T* myPsi_gs  = new Psi_T(s_min, m_min);
    double value;
    double value2;
    double error;
    int n_bins=500;
    double bins_size= double(xmax)/double(n_bins);
    vector<int> walker(n_bins);
    L=10000;
    
    for(int k=0; k<1000; k++){ //equilibration
        x = Metropolis_psi(x, myPsi_gs, rnd, nhit, ntot);
    }
    for(int n=0; n<N; n++){
        integral=0;
        nhit=0;
        ntot=0;
        for(int i=0; i<L; i++){
            x = Metropolis_psi(x, myPsi_gs, rnd, nhit, ntot);
            walker[int( (x+xmax)/(2.*bins_size))] += 1;
            integral += (-0.5*(myPsi_gs->Second_derivate(x)) + (myPot->Eval(x))*(myPsi_gs->Eval(x)) ) /(myPsi_gs->Eval(x));
        }
        integral=integral/double(L);
        value=(value*n+integral)/double(n+1);
        value2=(value2*n+integral*integral)/double(n+1);
        error= Error(value, value2, n);
        gs_eigenvalue << n << " " << value << " " << value2 << " " << error << " " << integral << endl;
    }
    for(int i=0; i< n_bins; i++){
        myfile << -xmax + xmax*i*2/double(n_bins) + xmax/double(2*n_bins) << " " << walker[i]/double(11000) << endl;
    }
    cout << "++++++++++  " << xmax/double(n_bins) << "++++++++++++" << endl;
    
    myfile.close();
    gs_eigenvalue.close();

 
    rnd->SaveSeed();
    return 0;
}




    
