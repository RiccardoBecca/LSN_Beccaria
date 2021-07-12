
#ifndef __Random__
#define __Random__
#include <vector>
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};

double Calcoloerrore(double, double);

double Error(std::vector<double>, std::vector<double>, int n);
/*
double error_becca (double x, double y, int n){
    if(n==0){
        return 0;
    }
    else{
        return sqrt((y- pow(x,2))/double(n));
    }
};

double error_becca(){
};
*/
/*
void func(std::vector<int> vect){
    if(vect.size()!=0){
        std::cout << vect[0] << std::endl;
    }

}
*/
#endif // __Random__
/*
class errore {

public:
    double Valore(vector<double> x, vector<double> y, int n);
};

*/
