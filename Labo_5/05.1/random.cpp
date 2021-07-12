

#include "random.h"


using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double lambda){
    double y=Rannyu();
    return -(log(1-y))/(lambda);
}

double Random :: Lorentz(double mu, double gamma){
    double y=Rannyu();
    return mu +gamma*tan(M_PI*(y-0.5));
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Pi(){
    double y = Rannyu();
    double x = Rannyu();
    if(x==0){
        return 0.;
    }
    if(sqrt(x*x+y*y)<1){
        return(2*atan(y/x));
    }
    else{
        return Pi();
    }
}


double Random :: d_prob(double a, double b, double f_max, FunzioneBase* f){
    double x = Rannyu(a,b);
    double r = Rannyu();
    if (r<((f->Eval(x))/f_max)){
        return x;
    }
    else {
        return d_prob(a, b, f_max, f);
    }
    
}


void Random :: Passo(vector<double> & v){
    if(v.size()!=3){
        cout << "There is something wrong mate!" << endl << "Quit program..." << endl;
        exit(0);
    }
    else{
        int i=Rannyu(0,3);
        int s= (int)Rannyu(0,2)*2-1;
        v[i]+=s;
    }
}

void Random :: Passo_continuo(vector<double> & v){
    if(v.size()!=3){
        cout << "There is something wrong mate!" << endl << "Quit program..." << endl;
        exit(0);
    }
    else{
        double theta=acos(1-2*Rannyu());
        double fi= Rannyu(0,2*M_PI);
        
        v[0]+=sin(theta)*cos(fi);
        v[1]+=sin(theta)*sin(fi);
        v[2]+=cos(theta);
    }
}


double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}



