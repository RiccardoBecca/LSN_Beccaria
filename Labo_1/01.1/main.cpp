

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <iterator>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
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
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



    int N=150;
    int M=1000000;
    int L=double(M)/double(N);

    ofstream myfile;
    

    // Primo punto
    
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;
    
    for(int i=0; i<N;i++){
        double sum=0;
        for(int j=0;j<L;j++){
            sum += rnd.Rannyu();
        }
        double x=sum/double(L);
        
        
        if(i==0){
            sum_prog.push_back(x);
            su2_prog.push_back(x*x);
        }
        if(i!=0){
            sum_prog.push_back((sum_prog[i-1]*i+x)/double(i+1));
            su2_prog.push_back((su2_prog[i-1]*i+x*x)/double(i+1));
        }
        err_prog.push_back(Error(sum_prog, su2_prog, i));
        

    }

    myfile.open ("random_events.txt");
    
    for(int i=0; i< sum_prog.size(); i++){

        myfile << sum_prog[i] << " " << su2_prog[i] << " "<< err_prog[i] <<  endl ;
    }
    myfile.close();

    
    // Secondo punto

    vector<double> sum_prog2;
    vector<double> su2_prog2;
    vector<double> err_prog2;
    
    for(int i=0; i<N;i++){
        double sum=0;
        for(int j=0;j<L;j++){
            sum += pow(rnd.Rannyu()-0.5,2);
        }
        double x=sum/double(L);
        
        
        if(i==0){
            sum_prog2.push_back(x);
            su2_prog2.push_back(x*x);
        }
        if(i!=0){
            sum_prog2.push_back((sum_prog2[i-1]*i+x)/double(i+1));
            su2_prog2.push_back((su2_prog2[i-1]*i+x*x)/double(i+1));
        }
        err_prog2.push_back(Error(sum_prog2, su2_prog2, i));
        

    }

    myfile.open ("random_events_2.txt");
    
    for(int i=0; i< sum_prog2.size(); i++){

        myfile << sum_prog2[i] << " " << su2_prog2[i] << " "<< err_prog2[i] <<  endl ;
    }
    myfile.close();
    
    
    
    
    //TERZO PUNTO
    myfile.open ("chi_2.txt");
    for(int f=0; f<100; f++){
        int M2 = 100;
        int n=10000;
        vector<vector<double>> n_i( M2 , vector<double> (2, 0));
        for(double i=0;i<M2;i++){
            n_i[i][0]=double(i)/double(M2);
        }
  
        for(int j=0;j<n;j++){
            double x = rnd.Rannyu(0,1);
            if(x>0.99){
                n_i[99][1]+=1;
            }
            for(int t=0;t<M2-1;t++){
                if(n_i[t][0]<x && x<n_i[t+1][0]){
                    n_i[t][1]+=1;
                }
         
            }
         
        }
        double chi_2 = 0;
        for(int i=0; i<M2; i++){
            chi_2+=pow((n_i[i][1]-(double(n)/double(M2))),2)/(double(n)/double(M2));
        }

        myfile << chi_2 << endl;
        
    }
    myfile.close();
    
   rnd.SaveSeed();
   return 0;
}


