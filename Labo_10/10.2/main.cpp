#include "lib.h"
#include "mpi.h"
#include <string>

using namespace std;
 
int main (int argc, char *argv[]){
    
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    
    Random rnd;
    //Random* rnd=new Random();
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for(int t=0;t<rank+1;t++){
            Primes >> p1 >> p2;
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
    
    
    int n_generazioni, square_circle, n_cities, n_figli, n_migr;
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
    ReadInput >> n_migr;
    
    ReadInput.close();
    
    int rank1, rank2;
    int best_1 [n_cities];
    for(int y=0;y<n_cities;y++){
        best_1[y]=0;
    }
    int best_2 [n_cities];
    for(int y=0;y<n_cities;y++){
        best_2[y]=0;
    }
    
    ofstream best_len;
    string best_text=to_string(rank);
    best_text.append(".Best_len.txt");
    
    ofstream medium_len;
    string medium_text=to_string(rank);
    medium_text.append(".Medium_len.txt");
    
    best_len.open(best_text);
    medium_len.open(medium_text);

    Cities cities(square_circle, n_cities,rnd);
    Generazione mygen(n_cities, n_figli, prob_cross, prob_m1, prob_m2, prob_m3, prob_m4,rnd);
    mygen.Measure(cities);
    
    for(int i=0; i<n_generazioni; i++){
        if(i%n_migr ==0){
            cout << "GENERATION " << i << " out of 1000000 : " << 100*double(i)/double(n_generazioni) << " %" << endl;
        }

        
        if(i==0){
            ofstream out_first;
            string first_trip=to_string(rank);
            first_trip.append(".First_trip.txt");
            out_first.open(first_trip);
            for(int j=0; j<mygen.Get_son(0).Size(); j++){
                out_first << mygen.Get_son(0).Get(j) << " " << cities.Get(PBC(mygen.Get_son(0).Get(j) -1 ,n_cities),0) << " " << cities.Get(PBC(mygen.Get_son(0).Get(j) -1 ,n_cities),1) << endl;
            }
            out_first.close();
        }
        
        mygen.Measure(cities);
        mygen.Sorting();
        best_len << i << " " << mygen.Get_Best_Len() << endl;
        medium_len << i << " " << mygen.Get_Mean_Len() << endl;
        
        if(i%50 == 0){
            int itag=1;
                        
            if(rank==0){
                rank1= int(rnd.Rannyu(0.,size));
                rank2= int(rnd.Rannyu(0.,size));
                while(rank1==rank2){
                    rank2= int(rnd.Rannyu(0.,size));
                }
            }
            MPI_Bcast(&rank1, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&rank2, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (rank==rank1){
               //Itinerario Best1 = BestFit(pop);
                for (int t=0; t<n_cities; t++){ best_1[t]=mygen.Get_son(0).Get(t);
                }
            }
            if (rank==rank2){
                //cout << endl<< endl << "BEFORE " << endl;
                //mygen.Get_son(0).Print();
                for (int t=0; t<n_cities; t++){ best_2[t]=mygen.Get_son(0).Get(t);
                }
            }
            
            if(rank==rank1){
                MPI_Send(&best_1[0], n_cities, MPI_INT, rank2, itag, MPI_COMM_WORLD);
                MPI_Recv(&best_2[0], n_cities, MPI_INT, rank2, itag, MPI_COMM_WORLD, &stat);
                for(int r=0;r<n_cities;r++){
                    mygen.Get_son(0).Set(r, best_2[r]);
                }

            }
            
            if(rank==rank2){
                //cout <<endl << endl << " rank " << rank << " " << best_1[23]<< " " << best_1[24] << endl;
                
                MPI_Send(&best_2[0], n_cities, MPI_INT, rank1, itag, MPI_COMM_WORLD);
                MPI_Recv(&best_1[0], n_cities, MPI_INT, rank1, itag, MPI_COMM_WORLD, &stat);
                for(int r=0; r<n_cities;r++){
                    mygen.Get_son(0).Set(r, best_1[r]);
                }
                //cout << "AFTER " << endl;
                //for(int r=0; r<n_cities;r++){
                //    cout << best_1[r] << " ";
                //}
                //cout << endl;
                //mygen.Get_son(0).Print();
                //cout << endl << endl;
                
            }
            if(rank==rank2){
                Figlio f1 (n_cities);
                //cout << endl << "BEFORE " << endl;
                //mygen.Get_son(0).Print();
                //f1.Print();
                for(int r=0; r<n_cities;r++){
                    f1.Set(r,best_1[r]);
                }
                //cout << "AFTER" << endl;
                //f1.Print();
                //mygen.Get_son(0).Print();
                mygen.Change(f1,0);
                //mygen.Get_son(0).Print();
            }
            if(rank==rank1){
                Figlio f2(n_cities);
                for(int r=0; r<n_cities;r++){
                    f2.Set(r,best_2[r]);
                }
                mygen.Change(f2,0);
            }
            
            mygen.Measure(cities);
            mygen.Sorting();
        }

        
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
    string filename=to_string(rank);
    filename.append(".Final_trip.txt");
    out.open(filename);
    for(int i=0; i<mygen.Get_son(0).Size(); i++){
        out << mygen.Get_son(0).Get(i) << " " << cities.Get(PBC(mygen.Get_son(0).Get(i) -1 ,n_cities),0) << " " << cities.Get(PBC(mygen.Get_son(0).Get(i) -1 ,n_cities),1) << endl;
    }
    out.close();
 
    rnd.SaveSeed();
    //rnd->SaveSeed();
    MPI_Finalize();
    return 0;
    
}

