#include "lib.h"


using namespace std;
 
int main (int argc, char *argv[]){
    
    
    Random rnd;
    //Random* rnd=new Random();
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for(int j=0; j<9;j++){ //con 5 fin'ora è il migliore. Con 9 a meno di uno è uguale. 33 bellissimo
            Primes >> p1 >> p2 ;
        }

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
    
    
    int n_generazioni, square_circle, n_cities;
    double prob_cross, prob_m1, prob_m2, prob_m3, prob_m4;
    
    ifstream ReadInput;
    ReadInput.open("input.dat");
    
    ReadInput >> n_generazioni;
    ReadInput >> square_circle;//0 if cities on a cirle, 1 if cities into a square
    ReadInput >> n_cities;
    int n_figli=1;
    //int accepted=0;
    double T=3;
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
    }
    else{
        best_len.open("Best_len_square.txt");
    }
    
    medium_len.open("Medium_len.txt");
    
    Cities cities(square_circle, n_cities,rnd);
    Generazione mygen(n_cities, n_figli, prob_cross, prob_m1, prob_m2, prob_m3, prob_m4,rnd, 3);
    mygen.Measure(cities);
    
    while(T>0.3){
        cout << " I AM AT TEMPERATURE:  " << T << endl;
                
        if(T==3){
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
        
        for(int j=0;j<500;j++){
            //mygen.Metro(cities);
            mygen.Metro_swap(cities);
            mygen.Metro_permut(cities);
            mygen.Metro_shift(cities);
            mygen.Metro_inversion(cities);
        }
        T=T*0.997;
        mygen.Set_T(T);
        
        //mygen.Sorting();
        best_len << double(1./T) << " " << mygen.Get_Best_Len() << endl;

        double percent=(mygen.Get_accepted())/double(4*5000);
        cout << "PERCENTUAL METROPOLIS ACCEPTED : " << percent << endl;
        mygen.Set_accepted(0);
    }
     
    best_len.close();
    medium_len.close();
    //mygen.Sorting();
    
    
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

    return 0;
    
}

