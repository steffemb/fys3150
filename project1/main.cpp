#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "time.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
    {
    ofstream myfile;
    int n = atoi(argv[1]);
    //char *filename = new char[10000];
    //printf(filename, "file_%d.dat", n);
    //myfile.open(filename);
    myfile.open("file_100000.dat");
    clock_t start, finish; //declare start and final time
    start = clock();
    double a = -1;
    double b = 2;
    double c = -1;
    vec temp = zeros(n+2);
    vec f = zeros(n+2);
    vec u = zeros(n+2);
    double dx = 1.0/(n+1);
    vec btemp = zeros(n+2);

    for(int i=0 ; i < n+1 ; i++) {
        f(i) = dx*dx*100.0*exp(-10*dx*i);
    }
    btemp(1) = b;
    u[1] = f[1];

    for(int i=2 ; i <= n ; i++) {
        temp[i] = c/btemp[i-1];
        btemp[i] = b-a*temp[i];
        u[i] = f[i] - u[i-1]*temp[i];
    }
    u[n] = u[n] / btemp[n];

    for(int i=n-1 ; i >= 1 ; i--) {
        u[i] = (u[i] - u[i+1]*a)/btemp[i];
    }

    finish = clock();
    double totaltime = ( (finish - start)/((double)CLOCKS_PER_SEC) );

    //cout << totaltime << endl;
    {
    clock_t start, finish; //declare start and final time
    start = clock();

    mat A = zeros<mat>(n,n);
    A.diag(0).fill(2);
    A.diag(1).fill(-1);
    A.diag(-1).fill(-1);

    mat L, U;
    lu(L, U, A);
    vec f_without_end_points = f.subvec(1, n);

    vec y = solve(L,f_without_end_points);
    vec x = solve(U,y);

    finish = clock();
    double totaltime_LU = ( (finish - start)/((double)CLOCKS_PER_SEC) );


    for (int i=0; i<=n+1; i++) { // writes file with 3 columns, data for u''(x), x, and totaltime
        myfile << u[i] << ' ' << i*dx << ' ' << x[i] << ' ' << totaltime << ' ' << totaltime_LU << endl;
    //myfile << '0' << ' ' << '0' << ' ' << '0' << ' ' << totaltime << ' ' << totaltime_LU << endl;
    }

    myfile.close();
    }

    return 0;
}
