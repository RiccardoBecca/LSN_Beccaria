
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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


    
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    myfile.open("exp_central_theorem.txt");
    myfile2.open("lor_central_theorem.txt");
    myfile3.open("dice_central_theorem.txt");
    
    int n=pow(10,4);
    
    vector<int> N;
    N.push_back(1);
    N.push_back(2);
    N.push_back(10);
    N.push_back(100);
    

    for(int i=0; i<n;i++){
        for(int t=0; t<N.size();t++){
            double S_n=0;
            double S_n2=0;
            double S_n3=0;
            for(int j=0; j< N[t]; j++){
                S_n+=rnd.Exp(1);
                S_n2+=rnd.Lorentz(0.,1.);
                S_n3+=rnd.Rannyu();
            }
            S_n=S_n/double(N[t]);
            S_n2=S_n2/double(N[t]);
            S_n3=S_n3/double(N[t]);
            
            myfile << S_n << " ";
            myfile2 << S_n2 << " ";
            myfile3 << S_n3 << " ";
        }
        myfile << endl;
        myfile2 << endl;
        myfile3 << endl;
    }
    myfile.close();
    myfile2.close();
    myfile3.close();
    
    //myfile.open("distributions.txt");
    //for(int i=0; i<10000;i++){
    //    myfile << rnd.Exp(1.) << " " << rnd.Lorentz(3., 2.) << endl;
    //}
    
    //myfile.close();

    rnd.SaveSeed();
    return 0;
}


