#include "numericaldiff.h"
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

NumericalDiff::NumericalDiff()
{
}

void NumericalDiff::Analytical(mat &u, double dt, double dx, double dy){
    int N = 100;
    int x_steps = 1/dx;
    int y_steps = 1/dy;
    for(int i = 0; i<x_steps; i++){
        double x = i*dx;
        for(int j=0; j<y_steps; j++){
            double y = j*dy;
            //cout << "x = " << x << " y = " << y << " pi = " << M_PI << endl;
            for(int n=1; n<N; n++){
                u(j,i) += (2/(M_PI*n))*(1-cos(n*M_PI))*sin(n*M_PI*y)*( cosh(n*M_PI*x)- cosh(n*M_PI)/sinh(n*M_PI) * sinh(n*M_PI*x));
            }
        }
    }
}

void NumericalDiff::Explicit(mat &u, double dt, double dx, double dy){

    int x_steps = 1/dx;
    int y_steps = 1/dy;
    double h = dx; //problematic if dx != dy
    double alpha =dt/(h*h);
    mat nextu = u;

    for(int i = 1; i<x_steps-1; i++){
        for(int j=1; j<y_steps-1; j++){
            nextu(i,j) = u(i,j) + alpha*(u(i+1,j) +u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
        }
    }
    u = nextu;
}


void NumericalDiff::ImplicitJacobi(mat &u, int length, double h, double dt){
    mat nextu;
    double alpha,b;
    alpha = dt/(h*h); b = 1 / (1 + 4*alpha);
    nextu = u;
    for (int i=1; i<length-1; i++){
        for (int j=1; j<length-1; j++){
            nextu(i,j) = b*(u(i,j) + alpha*(nextu(i+1,j) + nextu(i-1,j) + nextu(i,j+1) + nextu(i,j-1)));
        }
    }
    u = nextu;
}
