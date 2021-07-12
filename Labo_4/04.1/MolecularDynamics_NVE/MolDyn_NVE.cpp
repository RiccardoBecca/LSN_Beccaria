
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
    if(start_measure != 1) nblk =1;
    nstep=nstep/nblk;
    for(int iblk=1; iblk <= nblk; ++iblk){
        if(start_measure==1) Reset(iblk);
    
        for(int istep=1; istep <= nstep; ++istep){
            Move();           //Move particles with Verlet algorithm
            if(istep%iprint == 0) cout << "Number of time-steps: " << istep*iblk << endl;
            if(istep%10 == 0){
                Measure();     //Properties measurement
                if(start_measure==1) Accumulate();
                //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf += 1;
            }
        }
        if(start_measure==1) Averages(iblk);
    }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.gas"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
    ReadInput >> nblk;
    ReadInput >> restart;
    ReadInput >> start_measure;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
    if(start_measure!=1) cout << "Number of steps = " << nstep << endl << endl;
    if(start_measure==1) cout << "Number of blocks = " << nblk << "; Step in every blocks = " << nstep/nblk << endl;
    
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
    
    
//measurement of g(r)
    igofr = 4;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box/2.0)/(double)nbins;


    if(restart==1){
        ReadConf.open("config.final");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
        
        ReadConf.open("config_pre.final");
        for (int i=0; i<npart; ++i){
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        
        Move();
        
        double v2=0.;
        
        for(int i=0; i<npart;i++){
            vx[i] = (Pbc(x[i]-xold[i]))/double(delta);
            vy[i] = (Pbc(y[i]-yold[i]))/double(delta);
            vz[i] = (Pbc(z[i]-zold[i]))/double(delta);
            v2+=vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        v2=v2/double(npart);
        
        double fs=sqrt((3.*temp)/v2);
        
        cout << "***********************************************************************" << endl;
        cout << " I am restarting not only at an actual configuration but also at an old " << endl;
        cout << " Actual Temperature T'=" << v2/double(3) << endl;
        cout << " Rescaling to Temperature T=" << temp << endl;
        cout << "Rescaling factor fs=" << fs << endl;
        cout << "***********************************************************************" << endl;
        
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }

    
    
    if(restart != 1){

        //Read initial configuration
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
        
        
        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
            }
    }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("./gas_phase/output_epot_gas.dat",ios::app);
  Ekin.open("./gas_phase/output_ekin_gas.dat",ios::app);
  Temp.open("./gas_phase/output_temp_gas.dat",ios::app);
  Etot.open("./gas_phase/output_etot_gas.dat",ios::app);
    
    
  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
        
//update of the histogram of g(r)
        walker[2+int(2*nbins*dr/double(box))] += 2;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }
    
    for (int k=igofr; k<igofr+nbins; ++k){ walker[k]=3.0*walker[k]/(rho*100*4*M_PI*(pow((k+1)*box/(2*nbins),3)- pow(k*box/(2*nbins),3)));
    }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    
    if(start_measure != 1){
        ofstream Equi;
        Equi.open("./gas_phase/equilibration_gas.dat",ios::app);
        Equi << stima_temp << endl;
    }
    
    
    walker[iv] = stima_pot;
      walker[it] = stima_temp;
      walker[ik]=stima_kin;
      walker[ie]=stima_etot;

    return;
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
    
   ofstream Pot_ave, Kin_ave, Temp_ave, Etot_ave, Gofr, Gave;
    double stima_ave_pot=0, stima_kin_ave=0, stima_temp_ave=0, stima_etot_ave=0;
    
    cout << "=*=*=*=*=*=*=*=*=*=" << endl;
    cout << " Block number " << iblk << endl;
    cout << "=*=*=*=*=*=*=*=*=*=" << endl;
    
    
    Pot_ave.open("./gas_phase/ave_epot_gas.out",ios::app);
    stima_ave_pot = blk_av[iv]/double(blk_norm); //Energy
    glob_av[iv]  += stima_ave_pot;
    glob_av2[iv] += stima_ave_pot*stima_ave_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    Pot_ave << iblk <<  " " << stima_ave_pot << " " << glob_av[iv]/(double)iblk << " " << err_pot << endl;
    Pot_ave.close();
    
    Kin_ave.open("./gas_phase/ave_kin_gas.out",ios::app);
    stima_kin_ave=blk_av[ik]/double(blk_norm);
    glob_av[ik]  += stima_kin_ave;
    glob_av2[ik] += stima_kin_ave*stima_kin_ave;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    Kin_ave << iblk <<  " " << stima_kin_ave << " " << glob_av[ik]/(double)iblk << " " << err_kin << endl;
    Kin_ave.close();
    
    Temp_ave.open("./gas_phase/ave_temp_gas.out",ios::app);
    stima_temp_ave=blk_av[it]/double(blk_norm);
    glob_av[it]  += stima_temp_ave;
    glob_av2[it] += stima_temp_ave*stima_temp_ave;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    
//    cout <<glob_av[it] << " "<< err_temp << " " << iblk << " "<< glob_av2[it] << endl;
    Temp_ave << iblk <<  " " << stima_temp_ave << " " << glob_av[it]/(double)iblk << " " << err_temp << endl;
    Temp_ave.close();
    
    Etot_ave.open("./gas_phase/ave_etot_gas.out",ios::app);
    stima_etot_ave=blk_av[ie]/double(blk_norm);
    glob_av[ie]  += stima_etot_ave;
    glob_av2[ie] += stima_etot_ave*stima_etot_ave;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
    Etot_ave << iblk <<  " " << stima_etot_ave << " " << glob_av[ie]/(double)iblk << " " << err_etot << endl;
    Etot_ave.close();
    
    Gofr.open("./gas_phase/output.gofr.0",ios::app);
    Gave.open("./gas_phase/output.gave.0",ios::app);
    
    for (int k=igofr; k<igofr+nbins; ++k){
        g = blk_av[k]/blk_norm;
        glob_av[k] += g;
        glob_av2[k] += g*g;
        err_gdir=Error(glob_av[k], glob_av2[k], iblk);
        Gofr << iblk << " " << (2*k+1)*box/double(4*nbins) << " " << glob_av[k]/(double)iblk << " " << err_gdir << endl;
        if(iblk==nblk){
            Gave << iblk << " " << (2*k+1)*box/double(4*nbins) << " " << glob_av[k]/(double)iblk << " " << err_gdir << endl;
        }
    }
    
    

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
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
    ofstream WriteConf_pre;

  cout << "Print final configuration to file config.final " << endl << endl;
    cout << "Print previous final configuration to config_pre.final" << endl << endl;
  WriteConf.open("config.final");
    WriteConf_pre.open("config_pre.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
      
      WriteConf_pre<< xold[i]/box << "    " << yold[i]/box << "    " << zold[i]/box << endl;
      
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
