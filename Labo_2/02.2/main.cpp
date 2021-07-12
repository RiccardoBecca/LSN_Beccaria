
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

    
    
    
    

    
    int M=10000;
    int N=100;
    int L=M/N;
    int passi=100;
    
    
    cout << endl << "Now performing lattice RW" << endl << "-------------------------" << endl << endl;
    
    ofstream myfile;
    myfile.open("random_walk.txt");
    
    vector<vector<double>> m( 100 , vector<double> (100, 0));

    for(int j=0; j<N; j++){
        cout << "Block " << j+1 << endl;
        for(int i=0; i<L; i++){
            vector<double> v (3,0);
            for(int t=0; t<passi; t++){
                m[j][t]+= (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
                rnd->Passo(v);
            }
        }
        for(int w=0; w<passi; w++){
            m[j][w]= sqrt(m[j][w]/double(L));
        }
    }
    
    //Average on blocks
    
    for(int t=0; t<passi; t++){
        double rn=0;
        double rn2=0;
        double e=0;
        for(int j=0; j<N; j++){
            rn=rn+m[j][t];
            rn2=rn2+m[j][t]*m[j][t];
        }
        rn=rn/double(N);
        rn2=rn2/double(N);
        e=sqrt((rn2-rn*rn)/double(N-1));
        
        myfile << t << " " << rn << " " << rn2 << " " << e <<endl;
    }
    
    myfile.close();
    
    
    cout << endl << endl << endl;
    
    
    
    cout << endl << "Now performing continuum RW" << endl << "-------------------------" << endl << endl;
    
    
    myfile.open("random_walk2.txt");
    
    vector<vector<double>> mc( 100 , vector<double> (100, 0));

    for(int j=0; j<N; j++){
        cout << "Block " << j+1 << endl;
        for(int i=0; i<L; i++){
            vector<double> v (3,0);
            for(int t=0; t<passi; t++){
                mc[j][t]+= (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
                rnd->Passo_continuo(v);
            }
        }
        for(int w=0; w<passi; w++){
            mc[j][w]= sqrt(mc[j][w]/double(L));
        }
    }
    
    //Average on blocks
    
    for(int t=0; t<passi; t++){
        double rn=0;
        double rn2=0;
        double e=0;
        for(int j=0; j<N; j++){
            rn=rn+mc[j][t];
            rn2=rn2+mc[j][t]*mc[j][t];
        }
        rn=rn/double(N);
        rn2=rn2/double(N);
        e=sqrt((rn2-rn*rn)/double(N-1));
        
        myfile << t << " " << rn << " " << rn2 << " " << e <<endl;
    }
    
    myfile.close();
    
    rnd->SaveSeed();
    return 0;
}


