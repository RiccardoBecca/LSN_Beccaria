
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <math.h>
#include "lib.h"
#include "random.h"

using namespace std;


double Error(std::vector<double> ave, std::vector<double> av2, int n){
    if(n==0){
        return 0;
    }
    else{
        return sqrt((av2[n]-ave[n]*ave[n])/double(n));
    }
}

double Error(double ave, double av2, int n){
    return sqrt((av2-ave*ave)/double(n));
}

double max(double a, double b){
    if (a<b){
        return b;
    }
    else{
        return a;
    }
}


