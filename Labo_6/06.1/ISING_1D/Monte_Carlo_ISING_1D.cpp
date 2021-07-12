

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
#include <vector>

using namespace std;

int main()
{

    
    
  Input(0); //Inizialization
    
    vector<double> sum_prog;
    vector<double> su2_prog;
    vector<double> err_prog;
    double x=0;
    double x2=0;
    
    
    // Example of equilibration
    temp=1.5;
    beta=1./1.5;
    ofstream myfile_equi;
    if(metro==1){
        myfile_equi.open("Equilibration_Metro.txt");
    }
    else{
        myfile_equi.open("Equilibration_Gibbs.txt");
    }
    if(temp==1.5){
        for(int t=0; t<100;t++){
            x=0;
            x2=0;
            for(int k=0; k<10000; k++){
                Move(metro);
                for (int i=0; i<nspin; ++i)
                {
                   x += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
                }
            }
            x=x/((double)nspin*10000);
            x2=x*x;
            if(t==0){
                sum_prog.push_back(x);
                su2_prog.push_back(x2);
            }
            if(t!=0){
                sum_prog.push_back((sum_prog[t-1]*t+x)/double(t+1));
                su2_prog.push_back((su2_prog[t-1]*t+x2)/double(t+1));
            }
            err_prog.push_back(Error(sum_prog, su2_prog, t));
        }
    }
    
    for(int i=0; i< sum_prog.size(); i++){

        myfile_equi << sum_prog[i] << " " << su2_prog[i] << " "<< err_prog[i] <<  endl ;
    }
    myfile_equi.close();
    
    temp=0.5;
    while (temp <=2.5){
        Input(1); //forse devo toglierlo
        
//        Next Lines are for equilibration
        
        cout << "                        START EQUILIBRATION                " << endl;
        accepted=0;
        attempted=0;
        for(int t=0; t<100000;t++){
            Move(metro);
        }
        
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move(metro);
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk);   //Print results for current block
        }
        temp =temp+ 0.1;
    }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(int a)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

//  ReadInput >> temp; //forse da scoommentare
    
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
    
  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();
    
    if(a==0){
        if(metro==1){
            remove("./Equilibration_Metro.txt");
            if(h==0.02){
                remove("./Myfinalfile_magnetization_Metro.txt");
            }
            else{
                remove("./Myfinalfile_heat_Metro.txt");
                remove("./Myfinalfile_ene_Metro.txt");
                remove("./Myfinalfile_susceptibility_Metro.txt");
            }
        }
        else{
            remove("./Equilibration_Gibbs.txt");
            if(h==0.02){
                remove("./Myfinalfile_magnetization_Gibbs.txt");
            }
            else{
                remove("./Myfinalfile_heat_Gibbs.txt");
                remove("./Myfinalfile_ene_Gibbs.txt");
                remove("./Myfinalfile_susceptibility_Gibbs.txt");
            }
        }
    }




//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  iu2= 4; //Energi^2
 
  n_props = 5; //Number of observables

//initial configuration
    if(restart==0){
        for (int i=0; i<nspin; ++i)
        {
          if(rnd.Rannyu() >= 0.5) s[i] = 1;
          else s[i] = -1;
        }
    }
    if(restart==1){
        ifstream final_conf;
        final_conf.open("config.final");
        for(int i=0; i<nspin;i++){
            final_conf >> s[i];
        }
        final_conf.close();
        
    }
  //for (int i=0; i<nspin; ++i)
  //{
  //  if(rnd.Rannyu() >= 0.5) s[i] = 1;
  //  else s[i] = -1;
  //}
     /*
    for (int i=0; i<nspin; ++i)
    {
      s[i] = -1;
    }
      */
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        double a=Min(1.,exp(-beta*(Boltzmann(-s[o],o)-Boltzmann(s[o], o))));
        double r=rnd.Rannyu(0,1);
        if(r<a){
            s[o]=-s[o];
            accepted +=1;
        }
        else{
            s[o]=s[o];
        }
        attempted+=1;
    }
    else //Gibbs sampling
    {
        double a= 1.0/(1.0 + exp (-2.0*beta*(J*(s[Pbc(o+1)]+s[Pbc(o-1)])+h)));
        double r=rnd.Rannyu(0,1);
        if(r<a){
            s[o]=1;
        }
        else{
            s[o]=-1;
        }
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0, u2=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
      u2+=pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]),2);
      m+=s[i];

// INCLUDE YOUR CODE HERE
  }
//    cout << m << " Hello " << endl;
  walker[iu] = u;
    walker[im] = m;
    walker[iu2]=u2;
    walker[ix]=m*m;
    
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(h==0){
    
        Ene.open("output.ene.0",ios::app);
        stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
        glob_av[iu]  += stima_u;
        glob_av2[iu] += stima_u*stima_u;
        err_u=Error(glob_av[iu],glob_av2[iu],iblk);
        Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
        Ene.close();
    
        if (iblk==20){
            ofstream myfinalfile_ene;
            if(metro==1){
                myfinalfile_ene.open("Myfinalfile_ene_Metro.txt", ios::app);
            }
            else{
                myfinalfile_ene.open("Myfinalfile_ene_Gibbs.txt", ios::app);
            }
            //myfinalfile_ene.open("Myfinalfile_ene.txt", ios::app);
            myfinalfile_ene << temp << " " << glob_av[iu]/(double)iblk << " " << err_c/double(100) << endl;
            glob_av[iu] = 0.;
            glob_av2[iu]= 0.;
        }
    }
    

// INCLUDE YOUR CODE HERE
    
    if(h==0){
    
        Heat.open("output.heat.0", ios::app);
        stima_u = blk_av[iu]/blk_norm/(double)nspin;
        stima_u2 = blk_av[iu2]/blk_norm/(double)nspin;

        glob_av[ic] += beta*beta*(stima_u2-stima_u*stima_u);
        glob_av2[ic]+= pow(-beta*beta*(stima_u2-stima_u*stima_u),2);
        err_c= Error(glob_av[ic],glob_av2[ic],iblk);
        Heat << setw(wd) << iblk <<  setw(wd) << -beta*beta*(stima_u2-stima_u*stima_u)/double(nspin) << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
        Heat.close();
    
        cout << "----------------------------" << endl << endl;
    
        if (iblk==20){
            ofstream myfinalfile_heat;
            if(metro==1){
                myfinalfile_heat.open("Myfinalfile_heat_Metro.txt", ios::app);
            }
            else{
                myfinalfile_heat.open("Myfinalfile_heat_Gibbs.txt", ios::app);
            }
            myfinalfile_heat << temp << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
            glob_av[ic] = 0.;
            glob_av2[ic]= 0.;
        }
    }
    
    if(h!=0){
    
        Mag.open("output.magnetization.0", ios::app);
        stima_m = blk_av[im]/blk_norm/(double)nspin;
        glob_av[im] += stima_m;
        glob_av2[im]+= stima_m*stima_m;
        err_c= Error(glob_av[im],glob_av2[im],iblk);
        Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_c << endl;
        Mag.close();
    
        cout << "----------------------------" << endl << endl;

        if (iblk==20){
            ofstream myfinalfile_mag;
            if(metro==1){
                myfinalfile_mag.open("Myfinalfile_magnetization_Metro.txt", ios::app);
            }
            else{
                myfinalfile_mag.open("Myfinalfile_magnetization_Gibbs.txt", ios::app);
            }
            //myfinalfile_mag.open("Myfinalfile_magnetization.txt", ios::app);
            myfinalfile_mag << temp << " " << glob_av[im]/(double)iblk << " " << err_c << endl;
            glob_av[im] = 0.;
            glob_av2[im]= 0.;
        
        }
    }
    
    if(h==0){
    
        Chi.open("output.susceptibility.0", ios::app);
        stima_x = blk_av[ix]/blk_norm/(double)nspin;
        glob_av[ix] += beta*stima_x;
        glob_av2[ix]+= pow(beta*stima_x,2);
        err_c= Error(glob_av[ix],glob_av2[ix],iblk);
        Chi << setw(wd) << iblk <<  setw(wd) << beta*stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_c << endl;
        Chi.close();
    
        if (iblk==20){
            ofstream myfinalfile_chi;
            if(metro==1){
                myfinalfile_chi.open("Myfinalfile_susceptibility_Metro.txt", ios::app);
            }
            else{
                myfinalfile_chi.open("Myfinalfile_susceptibility_Gibbs.txt", ios::app);
            }
            //myfinalfile_chi.open("Myfinalfile_susceptibility.txt", ios::app);
            myfinalfile_chi << temp << " " << glob_av[ix]/(double)iblk << " " << err_c << endl;
            glob_av[ix] = 0.;
            glob_av2[ix]= 0.;
        
        }
    }
    
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


double Min(double a, double b)
{
    if (a<b){
        return a;
    }
    else{
        return b;
    }
}

double Error(vector<double> ave, vector<double> av2, int n){
    if(n==0){
        return 0;
    }
    else{
        return sqrt((av2[n]-ave[n]*ave[n])/double(n));
    }
}
