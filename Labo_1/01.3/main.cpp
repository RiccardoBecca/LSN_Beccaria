
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "lib.h"

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


    
    int M=3000000;
    int N=200;
    int L=M/N;
    
    double l=1.;
    double d=1.2;
    
    
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;


    
    ofstream myfile;
    myfile.open("pi.txt");
    
    for(int i=0; i<N;i++){
        double media=0.;
        for(int j=0; j<L; j++){
            double y=rnd.Rannyu();
            double theta=rnd.Pi();
            double y_max=y*d+l*sin(theta)/double(2);
            if(y_max>=d){
                media+=1;
            }
            double y_min=y*d-l*sin(theta)/double(2);
            if(y_min<=0){
                media+=1;
            }
        }
        media=media/double(L);
        double x=2*l/double(d*media);
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
    for(int i=0; i< sum_prog.size(); i++){

        myfile << sum_prog[i] << " " << su2_prog[i] << " " << err_prog[i] << endl ;
    }
    
    myfile.close();
 


    rnd.SaveSeed();
    return 0;
}


