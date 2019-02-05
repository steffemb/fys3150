//#include <iostream>
//#include <armadillo>
//#include <fstream>
//#include <cmath>
//#include <cstdlib>
//#include "time.h"

//using namespace std;
//using namespace arma;

//int main(int argc, char* argv[])
//    {
//    ofstream myfile;
//    int n = atoi(argv[1]);
//    myfile.open("file_10.dat");
//    clock_t start, finish; //declare start and final time
//    start = clock();

//    vec f = zeros(n+2);
//    double dx = 1.0/(n+1);

//    for(int i=0 ; i < n+1 ; i++) {
//        f(i) = 100.0*exp(-10*dx*i);
//    }

//    mat A = eye<mat>(n,n);
//    mat x = solve(A, f);

//    finish = clock();
//    double totaltime = ( (finish - start)/((double)CLOCKS_PER_SEC) );

//    //for (int i=0; i<=n+1; i++) {
//    //    myfile << u_LU[i] << endl;
//    //}

//    myfile.close();

//    return 0;
//}
