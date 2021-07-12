#include "lib.h"


using namespace std;
 
int main (int argc, char *argv[]){
    
    
    Random rnd;
    //Random* rnd=new Random();
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
             //rnd->SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    
    int n_generazioni, square_circle, n_cities, n_figli;
    double prob_cross, prob_m1, prob_m2, prob_m3, prob_m4;
    
    ifstream ReadInput;
    ReadInput.open("input.dat");
    
    ReadInput >> n_generazioni;
    ReadInput >> square_circle;//0 if cities on a cirle, 1 if cities into a square
    ReadInput >> n_cities;
    ReadInput >> n_figli;
    ReadInput >> prob_cross;
    ReadInput >> prob_m1;//0 if cities on a cirle, 1 if cities into a square
    ReadInput >> prob_m2;
    ReadInput >> prob_m3;
    ReadInput >> prob_m4;
    
    ReadInput.close();
    
    
    ofstream best_len;
    ofstream medium_len;
    ofstream out_cities;
    if(square_circle==0){
        best_len.open("Best_len_circle.txt");
        medium_len.open("Medium_len_circle.txt");
        out_cities.open("Cities_circle.txt");
    }
    else{
        best_len.open("Best_len_square.txt");
        medium_len.open("Medium_len_square.txt");
        out_cities.open("Cities_square.txt");
    }

    Cities cities(square_circle, n_cities,rnd);
    for(int k=0; k<n_cities;k++){
        for(int z=0; z<2;z++){
            out_cities << cities.Get(k,z) << " ";
        }
        out_cities << endl;
    }
    out_cities.close();
    Generazione mygen(n_cities, n_figli, prob_cross, prob_m1, prob_m2, prob_m3, prob_m4,rnd);
    mygen.Measure(cities);
    
    for(int i=0; i<n_generazioni; i++){
        if(i%10 ==0){
            cout << "GENERATION " << i << " out of 1000000 : " << 100*double(i)/double(n_generazioni) << " %" << endl;
        }

        
        if(i==0){
            ofstream out_first;
            if(square_circle==0){
                out_first.open("First_trip_circle.txt");
            }
            else{
                out_first.open("First_trip_square.txt");
            }
            for(int j=0; j<mygen.Get_son(0).Size(); j++){
                out_first << mygen.Get_son(0).Get(j) << " " << cities.Get(PBC(mygen.Get_son(0).Get(j) -1 ,n_cities),0) << " " << cities.Get(PBC(mygen.Get_son(0).Get(j) -1 ,n_cities),1) << endl;
            }
            out_first.close();
        }

        mygen.Sorting();
        best_len << i << " " << mygen.Get_Best_Len() << endl;
        medium_len << i << " " << mygen.Get_Mean_Len() << endl;
        
        Generazione oldgen(n_cities,n_figli, prob_cross, prob_m1, prob_m2, prob_m3, prob_m4,rnd);
        oldgen.Copy(mygen);
        oldgen.Crossing(cities, mygen);
        mygen.Copy(oldgen);
        mygen.Mutation(cities);
        
    }
     
    best_len.close();
    medium_len.close();
    mygen.Sorting();
    
    
    ofstream out;
    if(square_circle==0){
        out.open("Final_trip_circle.txt");
    }
    else{
        out.open("Final_trip_square.txt");
    }
    for(int i=0; i<mygen.Get_son(0).Size(); i++){
        out << mygen.Get_son(0).Get(i) << " " << cities.Get(PBC(mygen.Get_son(0).Get(i) -1 ,n_cities),0) << " " << cities.Get(PBC(mygen.Get_son(0).Get(i) -1 ,n_cities),1) << endl;
    }
    out.close();
 
    rnd.SaveSeed();
    //rnd->SaveSeed();
    return 0;
    
}

