#include <iostream>
#include <armadillo>
#include <solver.h>
#include <omp.h>

using namespace std;
using namespace arma;



int main()
{
    ofstream myfile;
    myfile.open("data_time.dat");
    double i = 10; // spatial points
    double j = 100; // timesteps
    double l0 = 0.1;
    double dt = (l0*l0)/(2);
    //double dt = 0.04;
    double dx = 1.0/(i-1);
    cout << "alpha = " <<  dt/(dx*dx) << " dt = " << dt << endl;

    //Forward Euler:
    vec initialx = zeros(i);//linspace(0,1,i);
    initialx(0) = 1; //all particles at x = 0
    mat u = zeros<mat>(i,j);
    for(int k = 0; k<j; k++){
        u.col(k) = initialx;
    }

    //cout << u << endl;


    double wtime;
    wtime = omp_get_wtime();
    Solver mySolver; //choose method:
    //mySolver.ForwardEuler(u, i, j, dt, dx); //Explicit
    //mySolver.BackwardEuler(u, i, j, dt, dx); //Implicit
    mySolver.CrankNichols(u, i, j, dt, dx); //CranckNichols scheme
    //cout << u << endl;
    wtime = omp_get_wtime() - wtime;
    cout << "Elapsed seconds = " << wtime << endl;

    for(int  k = 0; k < i; k++){
        myfile << u.row(k);
    }
    //myfile << u << endl;
    myfile.close();
    //cout << u << endl;

    return 0;
}

