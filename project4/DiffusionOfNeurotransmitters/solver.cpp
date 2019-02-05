#include "solver.h"
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;

Solver::Solver()
{
}

void Solver::ClosedForm(mat &u, double n, double m, double dt, double dx){ //closed form solution with L = 1, and n = 1. (fys 3150 see lecture notes page 314, 2014)
    //vec x = linspace(0,1,n);

    for(int  j = 0; j < m; j++){ //needs initial values at u0
        double t = j*dt;
        for(int i = 0; i < n; i++){
            double x = i*dx;
            double a = M_PI*x;
            double A = 2*( (sin(a*x) - a*(x-1)*cos(a*x) )/(a*a) ); //-1/sqrt(4*M_PI);
            u(i,j) = -(1-(1-exp(-t))*x-exp(-x*t))+(1-x);//A*sin(a)*exp(-M_PI*M_PI*t) + (1-x); //(1/( sqrt(4*M_PI)) )*exp(- ((x(i)*x(i))/4*t) );
        }
    }
}


void Solver::ForwardEuler(mat &u, double n, double m, double dt, double dx){
    for(int  j = 0; j < m-1; j++){ //needs initial values at u0
        for(int i = 1; i < n-1; i++){
            u(i,j+1) = ( (u(i+1,j) - 2*u(i,j) + u(i-1,j))*(dt/(dx*dx)) + u(i,j) );
        }

    }

}

void Solver::BackwardEuler(mat &u, double n, double m, double dt, double dx){

    double alpha = dt/(dx*dx);
    double a = -alpha;
    double b = 1.0+2.0*alpha;
    double c = -alpha;
    cout << u.col(0);
    for(int  j = 1; j < m-1; j++){
        //vec temp = zeros(n);
        vec v = u.col(j-1);
        v(1)+=alpha;
        vec btemp = zeros(n);

        btemp(1) = b;
        //v[0] = 0;


        vec p =zeros<vec>(n);
        //forward substutution
        for(int i=2 ; i < n-1 ; i++) {
            double temp = c/btemp[i-1];
            btemp[i] = b-a*temp;
            v[i] = v[i] - v[i-1]*temp;
            cout << "forste lokke" << endl;
            cout << v[i] << endl;
            cout << btemp[i] << endl;
        }
        p[n-2] = v[n-2] / btemp[n-2];

        //backward substitution
        for(int i=n-3 ; i > 0 ; i--) {
            p[i] = (v[i] - p[i+1]*a)/btemp[i];
            cout << "andre lokke" << endl;
            cout << v[i] << endl;

        }
        //cout << v << endl;
        p(0) = 1;
        u.col(j) = p;
    }
}

void Solver::CrankNichols(mat &u, double n, double m, double dt, double dx){
    double alpha = dt/(dx*dx);
    double a = -alpha;
    double b = 2.0+2.0*alpha;
    double c = -alpha;
    cout << u.col(0);

    for(int  j = 1; j < m-1; j++){
        vec vEul = zeros(n);
        for(int i = 1; i < n-1; i++){
            vEul(i) = ( (u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))*(alpha) + 2*u(i,j-1) );
        }
        //vec temp = zeros(n);
        vec vImpl = vEul; //u.col(j-1);
        vImpl(1)+=alpha;
        vec btemp = zeros(n);

        btemp(1) = b;
        //v[0] = 0;


        vec p =zeros<vec>(n);
        //forward substutution
        for(int i=2 ; i < n-1 ; i++) {
            double temp = c/btemp[i-1];
            btemp[i] = b-a*temp;
            vImpl[i] = vImpl[i] - vImpl[i-1]*temp;
            cout << "forste lokke" << endl;
            cout << vImpl[i] << endl;
            cout << btemp[i] << endl;
        }
        p[n-2] = vImpl[n-2] / btemp[n-2];

        //backward substitution
        for(int i=n-3 ; i > 0 ; i--) {
            p[i] = (vImpl[i] - p[i+1]*a)/btemp[i];
            cout << "andre lokke" << endl;
            cout << vImpl[i] << endl;

        }
        //cout << v << endl;
        p(0) = 1;
        u.col(j) = (p);
    }
}

