
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

    
    
    int M=1000000;
    int N=1000;
    int L=M/N;
    
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;
    
    vector<double> sum_put_prog;
    vector<double> su2_put_prog;
    vector<double> err_put_prog;
    
    double S_in=100.;
    double T=1.;
    double K=100.;
    double r=0.1;
    double s=0.25;
    
    double S_i=0.;
    double z_i=0.;
    double m_call=0.;
    double m_put=0.;
    
    ofstream myfile_call;
    myfile_call.open("Call_cost.txt");
    
    ofstream myfile_put;
    myfile_put.open("Put_cost.txt");
    cout << " I am sampling directly the final asset price S(T)! " << endl;
    
    for(int j=0; j<N; j++){
        double Cn_call=0;
        double Cn_put=0;
        for(int i=0; i<L; i++){
            z_i= rnd->Gauss(0,1);
            S_i= S_in*exp((r-s*s*0.5)*T+s*z_i*sqrt(T));
            m_call=max(0.,S_i-K);
            m_put= max(0.,K-S_i);
            
            Cn_call+=exp(-r*T)*m_call;
            Cn_put+=exp(-r*T)*m_put;
        }
        Cn_call=Cn_call/double(L);
        Cn_put=Cn_put/double(L);
        if(j==0){
            sum_prog.push_back(Cn_call);
            su2_prog.push_back(Cn_call*Cn_call);
            
            sum_put_prog.push_back(Cn_put);
            su2_put_prog.push_back(Cn_put*Cn_put);
        }
        else{
            sum_prog.push_back((sum_prog[j-1]*j+Cn_call)/double(j+1));
            su2_prog.push_back((su2_prog[j-1]*j+Cn_call*Cn_call)/double(j+1));
            
            sum_put_prog.push_back((sum_put_prog[j-1]*j+Cn_put)/double(j+1));
            su2_put_prog.push_back((su2_put_prog[j-1]*j+Cn_put*Cn_put)/double(j+1));
            
        }
        err_prog.push_back(Error(sum_prog, su2_prog,j));
        err_put_prog.push_back(Error(sum_put_prog, su2_put_prog,j));
    }

    
    for(int i=0; i< sum_prog.size(); i++){
        myfile_call << sum_prog[i] << " " << su2_prog[i] << " " << err_prog[i] << endl;
        myfile_put << sum_put_prog[i] << " " << su2_put_prog[i] << " " << err_put_prog[i] << endl;
    }
    
    
    myfile_call.close();
    myfile_put.close();
    
    
    vector<double> Sum_prog;
    vector<double> Su2_prog;
    vector<double> Err_prog;
    
    vector<double> Sum_put_prog;
    vector<double> Su2_put_prog;
    vector<double> Err_put_prog;
    
    double S_ti=0.;
    double S_ti1=0.;
    double W=100;
    
    
    
    myfile_call.open("Call_cost_discrete.txt");
    myfile_put.open("Put_cost_discrete.txt");
    
    cout << " I am sampling the discretized path of the final asset price S(T)! " << endl;
    
    for(int j=0; j<N; j++){
        double Cn_call=0;
        double Cn_put=0;
        for(int i=0; i<L; i++){
            S_ti1=S_in;
            for(int w=0; w<W; w++){
                S_ti=S_ti1;
                z_i= rnd->Gauss(0,1);
                S_ti1= S_ti*exp((r-s*s*0.5)*(0.01)+s*z_i*sqrt(0.01));
            }
            m_call=max(0.,S_ti1-K);
            m_put= max(0.,K-S_ti1);
            Cn_call+=exp(-r*T)*m_call;
            Cn_put+=exp(-r*T)*m_put;
        }
        
        
        
        Cn_call=Cn_call/double(L);
        Cn_put=Cn_put/double(L);
        
        if(j==0){
            Sum_prog.push_back(Cn_call);
            Su2_prog.push_back(Cn_call*Cn_call);
            
            Sum_put_prog.push_back(Cn_put);
            Su2_put_prog.push_back(Cn_put*Cn_put);
        }
        else{
            Sum_prog.push_back((Sum_prog[j-1]*j+Cn_call)/double(j+1));
            Su2_prog.push_back((Su2_prog[j-1]*j+Cn_call*Cn_call)/double(j+1));
            
            Sum_put_prog.push_back((Sum_put_prog[j-1]*j+Cn_put)/double(j+1));
            Su2_put_prog.push_back((Su2_put_prog[j-1]*j+Cn_put*Cn_put)/double(j+1));
        }
        Err_prog.push_back(Error(Sum_prog, Su2_prog,j));
        Err_put_prog.push_back(Error(Sum_put_prog, Su2_put_prog,j));
    }

    
    for(int i=0; i< Sum_prog.size(); i++){
        myfile_call << Sum_prog[i] << " " << Su2_prog[i] << " " << Err_prog[i] << endl;
        myfile_put << Sum_put_prog[i] << " " << Su2_put_prog[i] << " " << Err_put_prog[i] << endl;
    }
    
    
    myfile_call.close();
    myfile_call.close();
    
    
    rnd->SaveSeed();
    return 0;
}




    
