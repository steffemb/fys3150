#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <cstdlib>

//#include "time.h"

using namespace std;
using namespace arma;

void findMaxValueInMatrix(mat &a, double &max_a, int &index_k, int &index_l) {
    int n = a.n_cols;
    max_a = 0;
    for(int k = 0 ; k < n ; k++) {// find a_kl
        for(int l = k + 1; l < n ; l++) {
            if(a(k,l)*a(k,l) > max_a){
                max_a = a(k,l)*a(k,l);
                index_k = k;
                index_l = l;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    ofstream myfile;
    int n = atoi(argv[1]);
    //char *filename = new char[10000]; //can i get lager n?
    //printf(filename, "file_%d.dat", n);
    //myfile.open(filename);
    myfile.open("eigenvec.dat");
    double epsilon = pow(10,-8);
    int index_l = -1;
    int index_k = -1;
    double max_a = 0;
    mat a = zeros<mat>(n, n);
    mat b = zeros<mat>(n, n);
    double tau;
    double t;
    vec ro = zeros(n);
    double ro_min = 0;
    double ro_max = 5.5;
    double h = (ro_max-ro_min)/(n+1);
    double e = -(1/(h*h));
    vec V = zeros(n);
    vec d = zeros(n);
    double r_ik;
    double r_il;
    mat R = zeros<mat>(n, n);
    double LowestLamda = 100000;
    double LowestLamdaIndex = 0;
    double NextLowest, ThirdLowest;
    double omega = 0.01;

    //fill in matrix a (schrödinger eq)
    a.diag(1).fill(e);
    a.diag(-1).fill(e);
    for(int i = 0; i < n; i++) {
        ro(i) = ro_min + (i+1)*h;
        //V(i) = ro(i)*ro(i); // for regular ocillator potential
        V(i) = omega*omega*ro(i)*ro(i) + 1/ro(i);
        d(i) = 2/(h*h) + V(i);
        a(i,i) = d(i);
    }

    // Setting up the eigenvector matrix
    for(int i = 0; i < n; i++ ){
        for(int j = 0; j < n; j++ ){
            if( i == j ){
                R(i, j) = 1.0;
            }else{
                R(i, j) = 0.0;
            }
        }
    }


    vec eigval;  // armadillo real eigenvalue/eigen vector
    mat eigvec;
    eig_sym( eigval, eigvec, a );
    cout << eigval.subvec(0,2) << endl;
    //cout << eigvec << endl;

    b=a;
    findMaxValueInMatrix(a, max_a, index_k, index_l);
    while(max_a > epsilon) {
        tau = (a(index_l,index_l) - a(index_k,index_k))/(2*a(index_k,index_l));
        if( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }else{
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }


        double c = 1/sqrt(1 + t*t); // cos(tetha)
        double s = t*c; // sin(theta)

        for(int i = 0 ; i < n ; i++) {// find a_kl
            if(i != index_k){
                if(i != index_l){
                    b(i,index_k) = a(i,index_k)*c - a(i,index_l)*s;
                    b(i,index_l) = a(i,index_l)*c + a(i,index_k)*s;
                    b(index_k,i) = a(i,index_k)*c - a(i,index_l)*s;
                    b(index_l,i) = a(i,index_l)*c + a(i,index_k)*s;
                }
            }
        }
        b(index_k,index_k) = a(index_k,index_k)*c*c - 2*a(index_k,index_l)*c*s + a(index_l,index_l)*s*s;
        b(index_l,index_l) = a(index_l,index_l)*c*c + 2*a(index_k,index_l)*c*s + a(index_k,index_k)*s*s;
        b(index_k,index_l) = 0;//(a(index_k,index_k) - a(index_l,index_l))*c*s + a(index_k,index_l)*(c*c - s*s);
        b(index_l,index_k) = 0;//(a(index_k,index_k) - a(index_l,index_l))*c*s + a(index_k,index_l)*(c*c - s*s);

        a = b;
        // compute the new eigenvectors
        for(int i = 0 ; i < n ; i++) {
            r_ik = R(i, index_k);
            r_il = R(i, index_l);
            R(i, index_k) = c*r_ik - s*r_il;
            R(i, index_l) = c*r_il + s*r_ik;
        }

        findMaxValueInMatrix(a, max_a, index_k, index_l);
        max_a = sqrt(max_a);



    }

    // writes file with 3 columns, data for u''(x), x, and totaltime
    //myfile << b(i,i) << endl; //' ' << i*dx << ' ' << x[i] << ' ' << totaltime << ' ' << totaltime_LU << endl;
    vec b_vec = b.diag();
    umat LamdaIndex = sort_index(b_vec);

    LowestLamda = LamdaIndex(0);
    NextLowest = LamdaIndex(1);
    ThirdLowest = LamdaIndex(2);
    //cout << b_vec(0) << b_vec(1) << b_vec(2) << endl;

    b_vec = sort(b_vec);
    cout << b_vec.subvec(0,2) << endl;

    //cout << LowestLamda << endl;
    for (int i=0; i < n; i++) {
        myfile << R(i, LowestLamda) << ' ' << R(i, NextLowest) << ' ' << R(i, ThirdLowest) << ' ' << ro_max << endl;
    }
    //cout << R << endl;
    myfile.close();

    return 0;
}

