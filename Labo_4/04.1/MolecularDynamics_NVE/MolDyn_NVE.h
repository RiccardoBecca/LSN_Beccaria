
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;
double walker[m_props];

double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_ave_pot;
double err_pot, err_kin, err_temp, err_etot, err_gdir;

//g(r)
int igofr;
double bin_size,nbins,sd,g;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];


// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, restart, nblk, start_measure;
double delta;

//functions
void Input(void);
void Move(void);
void Accumulate(void);
void Averages(int);
void Reset(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
