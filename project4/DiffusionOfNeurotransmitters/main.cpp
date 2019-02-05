#include <iostream>
#include <armadillo>
#include <solver.h>

using namespace std;
using namespace arma;



int main()
{
    ofstream myfile;
    myfile.open("data_time.dat");
    double i = 10; // spatial points
    double j = 100; // timesteps
    double dt = 0.04;
    double dx = 1.0/(i-1);
    cout << dt/(dx*dx) << endl;

    //Forward Euler:
    vec initialx = zeros(i);//linspace(0,1,i);
    initialx(0) = 1; //all particles at x = 0
    mat u = zeros<mat>(i,j);
    for(int k = 0; k<j; k++){
        u.col(k) = initialx;
    }

    //cout << u << endl;


    Solver mySolver; //choose method:
    //mySolver.ForwardEuler(u, i, j, dt, dx); //Explicit
    //mySolver.BackwardEuler(u, i, j, dt, dx); //Implicit
    //mySolver.CrankNichols(u, i, j, dt, dx); //CranckNichols scheme
    mySolver.ClosedForm(u, i, j, dt, dx);
    //cout << u << endl;
    for(int  k = 0; k < i; k++){
        myfile << u.row(k);
    }
    //myfile << u << endl;
    myfile.close();
    //cout << u << endl;

    return 0;
}

