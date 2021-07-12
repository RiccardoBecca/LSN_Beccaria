
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
    
    int n_hit=0;
    int n_tot=0;

  
    ofstream myfile;
    myfile.open("Medium_radius.txt");
    
    ofstream coordsfile;
    coordsfile.open("Coordinates_file.txt");
    
    
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;
    
    
    vector<double> v (3,0);
    H_gs_density* myH_gs  = new H_gs_density();
    double a=0.0529*pow(10,-9);
    double r=0;
    double s=a*1.2;
    
    v[0]=a;
    v[1]=a;
    v[2]=a;
    
    for(int i=0;i<N;i++){
        r=0;
        for(int j=0;j<L;j++){
            Metropolis_unif(v, myH_gs, rnd, s, n_hit, n_tot);
            r+=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            coordsfile << v[0] << " " << v[1] << " " << v[2] << endl;
        }
        r=r/double(L);
        if(i==0){
            sum_prog.push_back(r);
            su2_prog.push_back(r*r);
        }
        else{
            sum_prog.push_back((sum_prog[i-1]*i+r)/double(i+1));
            su2_prog.push_back((su2_prog[i-1]*i+r*r)/double(i+1));
        }
        err_prog.push_back(Error(sum_prog, su2_prog,i));
        myfile << sum_prog[i] << " " << su2_prog[i] << " "<< err_prog[i] << endl;
    }
    cout << "=========================================================" << endl;
    cout << "I'm using Metropolis with a uniform tentative probability" << endl << "with sigma = " << s << "  there's " << double(n_hit)*100/double(n_tot) <<"% of acceptance" << endl;
    cout << "(n, l, m) = (1, 0, 0)" << endl;
    cout << "=========================================================" << endl << endl << endl;

    
    myfile.close();
    coordsfile.close();
    
    ofstream my_gauss;
    my_gauss.open("Medium_radius_gauss.txt");
    
    n_hit=0;
    n_tot=0;
    
    vector<double> sum_gauss_prog;
    vector<double> su2_gauss_prog;
    vector<double> err_gauss_prog;
    
    
    vector<double> v2 (3,0);

    
    v2[0]=a;
    v2[1]=a;
    v2[2]=a;
    
    for(int i=0;i<N;i++){
        r=0;
        for(int j=0;j<L;j++){
            Metropolis_gauss(v2, myH_gs, rnd, s, n_hit, n_tot);
            r+=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
        }
        r=r/double(L);
        if(i==0){
            sum_gauss_prog.push_back(r);
            su2_gauss_prog.push_back(r*r);
        }
        else{
            sum_gauss_prog.push_back((sum_gauss_prog[i-1]*i+r)/double(i+1));
            su2_gauss_prog.push_back((su2_gauss_prog[i-1]*i+r*r)/double(i+1));
        }
        err_gauss_prog.push_back(Error(sum_gauss_prog, su2_gauss_prog,i));
        my_gauss << sum_gauss_prog[i] << " " << su2_gauss_prog[i] << " "<< err_gauss_prog[i] << endl;
    }
    cout << "==========================================================" << endl;
    cout << "I'm using Metropolis with a gaussian tentative probability" << endl << "with sigma = " << s << "  there's " << double(n_hit)*100/double(n_tot) <<"% of acceptance" << endl;
    cout << "(n, l, m) = (1, 0, 0)" << endl;
    cout << "==========================================================" << endl << endl << endl;
    my_gauss.close();
    
    
    
    // SECOND EIGENFUNCTION
     
    M=1000000;
    N=1000;
    L=M/N;
    
    n_hit=0;
    n_tot=0;

  
    myfile.open("Medium_radius_2.txt");
    
    coordsfile.open("Coordinates_file_2.txt");
    
    
    vector<double> sum_prog_2;
    vector<double> su2_prog_2;
    vector<double> err_prog_2;
    
    
    vector<double> v_2 (3,0);
    H_gs_density_2* myH_gs_2  = new H_gs_density_2();
    a=0.0529*pow(10,-9);
    r=0;
    s=a*1.2;
    
    v_2[0]=a;
    v_2[1]=a;
    v_2[2]=a;
    
    for(int i=0;i<N;i++){
        r=0;
        for(int j=0;j<L;j++){
            Metropolis_unif(v_2, myH_gs_2, rnd, s, n_hit, n_tot);
            r+=sqrt(v_2[0]*v_2[0]+v_2[1]*v_2[1]+v_2[2]*v_2[2]);
            //if(j%100==0){
            //    coordsfile << v_2[0] << " " << v_2[1] << " " << v_2[2] << endl;
            //}
            
        }
        r=r/double(L);
        if(i==0){
            sum_prog_2.push_back(r);
            su2_prog_2.push_back(r*r);
        }
        else{
            sum_prog_2.push_back((sum_prog_2[i-1]*i+r)/double(i+1));
            su2_prog_2.push_back((su2_prog_2[i-1]*i+r*r)/double(i+1));
        }
        err_prog_2.push_back(Error(sum_prog_2, su2_prog_2,i));
        myfile << sum_prog_2[i] << " " << su2_prog_2[i] << " "<< err_prog_2[i] << endl;
    }
    cout << "=========================================================" << endl;
    cout << "I'm using Metropolis with a uniform tentative probability" << endl << "with sigma = " << s << "  there's " << double(n_hit)*100/double(n_tot) <<"% of acceptance" << endl;
    cout << "(n, l, m) = (2, 1, 0)" << endl;
    cout << "=========================================================" << endl << endl << endl;

    
    myfile.close();
    //coordsfile.close();
    
    my_gauss.open("Medium_radius_gauss_2.txt");
    
    n_hit=0;
    n_tot=0;
    
    vector<double> sum_gauss_prog_2;
    vector<double> su2_gauss_prog_2;
    vector<double> err_gauss_prog_2;
    
    
    vector<double> v2_2 (3,0);

    
    v2_2[0]=a;
    v2_2[1]=a;
    v2_2[2]=a;
    
    for(int i=0;i<N;i++){
        r=0;
        for(int j=0;j<L;j++){
            Metropolis_gauss(v2_2, myH_gs_2, rnd, s, n_hit, n_tot);
            r+=sqrt(v2_2[0]*v2_2[0]+v2_2[1]*v2_2[1]+v2_2[2]*v2_2[2]);
            if(j%10==0){
                coordsfile << v2_2[0] << " " << v2_2[1] << " " << v2_2[2] << endl;
            }
        }
        r=r/double(L);
        if(i==0){
            sum_gauss_prog_2.push_back(r);
            su2_gauss_prog_2.push_back(r*r);
        }
        else{
            sum_gauss_prog_2.push_back((sum_gauss_prog_2[i-1]*i+r)/double(i+1));
            su2_gauss_prog_2.push_back((su2_gauss_prog_2[i-1]*i+r*r)/double(i+1));
        }
        err_gauss_prog_2.push_back(Error(sum_gauss_prog_2, su2_gauss_prog_2,i));
        my_gauss << sum_gauss_prog_2[i] << " " << su2_gauss_prog_2[i] << " "<< err_gauss_prog_2[i] << endl;
    }
    cout << "==========================================================" << endl;
    cout << "I'm using Metropolis with a gaussian tentative probability" << endl << "with sigma = " << s << "  there's " << double(n_hit)*100/double(n_tot) <<"% of acceptance" << endl;
    cout << "(n, l, m) = (2, 1, 0)" << endl;
    cout << "==========================================================" << endl << endl << endl;
    my_gauss.close();
    coordsfile.close();
    
    //SE SONO LONTANO
    
    n_hit=0;
    n_tot=0;

  
    myfile.open("Medium_radius_far.txt");
    
    coordsfile.open("Coordinates_file_far.txt");
    
    
    vector<double> sum_prog_far;
    vector<double> su2_prog_far;
    vector<double> err_prog_far;
    
    
    vector<double> v_far (3,0);
    H_gs_density* myH_gs_far  = new H_gs_density();
    a=0.0529*pow(10,-9);
    r=0;
    s=a*1.2;
    M=10000;
    N=100;
    L=M/N;
    
    v_far[0]=215*a;
    v_far[1]=215*a;
    v_far[2]=215*a;
    
    for(int i=0;i<N;i++){
        r=0;
        for(int j=0;j<L;j++){
            Metropolis_unif(v_far, myH_gs_far, rnd, s, n_hit, n_tot);
            r+=sqrt(v_far[0]*v_far[0]+v_far[1]*v_far[1]+v_far[2]*v_far[2]);
            coordsfile << v_far[0] << " " << v_far[1] << " " << v_far[2] << endl;
        }
        r=r/double(L);
        if(i==0){
            sum_prog_far.push_back(r);
            su2_prog_far.push_back(r*r);
        }
        else{
            sum_prog_far.push_back((sum_prog_far[i-1]*i+r)/double(i+1));
            su2_prog_far.push_back((su2_prog_far[i-1]*i+r*r)/double(i+1));
        }
        err_prog_far.push_back(Error(sum_prog_far, su2_prog_far,i));
        myfile << sum_prog_far[i] << " " << su2_prog_far[i] << " "<< err_prog_far[i] << endl;
    }
    cout << "=========================================================" << endl;
    cout << "STARTING POINT FAR FROM ORIGIN" << endl << "with sigma = " << s << "  there's " << double(n_hit)*100/double(n_tot) <<"% of acceptance" << endl;
    cout << "(n, l, m) = (1, 0, 0)" << endl;
    cout << "=========================================================" << endl << endl << endl;

    
    myfile.close();
    coordsfile.close();
    
    rnd->SaveSeed();
    return 0;
}




    
