#include <iostream>
#include <armadillo>
#include <vector>
#include <omp.h>
#include "particle.h"
#include "system.h"

#include <cmath>
#include <montecarlostep.h>
#include <meanvalue.h>
#include <numericaldiff.h>

using namespace arma;
using namespace std;




int main()
{

    /* -------------------------------------------------
     * Solving Diffusion equation in 2+1 dimensions:
     *
     --------------------------------------------------*/

    System mySystem1, mySystem2;
    MonteCarloStep myMonteCarlo;
    MeanValue myMeanValue;
    NumericalDiff myDiffSolver1, myDiffSolver2, myDiffSolver3;

    double N = pow(10, 4);                      //number of particles in system for MC methods
    double l0 = 0.1;                            // step length, dx = dy = l0
    double dt = (l0*l0)/(2);                    // invrted from formula for l0 for MC methods
    double T = 1.0;                             //final sim time
    int TimeSteps = T/dt;
    int length = int((1/l0))+1;
    double wtime1, wtime2, wtime3, wtime4;




    mat meanMat1 = zeros<mat>(int((1/l0)), int((1/l0)));

    for(int a = 0; a<int((1/l0)); a++){ //initial state for diff solver
        meanMat1(a,0) = 1;
    }

    mat meanMat5 = meanMat1;
    mat meanMat2 = meanMat1;                    // this is not necesarry for MC but we do it like this for simplicity.
    mat meanMat4 = meanMat1;

    cout << " timesteps = " << TimeSteps << endl;

//myMonteCarlo.RandNumTest();                   //for testing randomnumber generator
//    for(int i=0; i<N; i++){
//        //adding N particles to system
//        mySystem.addParticle(0.5);
//    }

    for(int k = 1; k<(int(1/l0))-1; k++){       //adding particles to system (setting up initial state)
        for(int i = 0; i<N/int(1/l0); i++){
            mySystem1.addParticle(l0*k+(l0/2));
            mySystem2.addParticle(l0*k+(l0/2));

        }
    }

    //opening files
    ofstream myfile1, myfile2, myfile3, myfile4;
    myfile1.open("MC_2d_gauss.dat"); myfile2.open("Explicit_diff_2d.dat"); myfile3.open("Analytical_diff_2d.dat");
    myfile4.open("MC_2d_fixed.dat");


    //one time loop for each method to compare the time each method uses
    //mark out the whole loop if the method is not needed
    //------------------------------------------------
    wtime1 = omp_get_wtime();
    for(int i = 0; i<TimeSteps; i++ ){
        myMonteCarlo.TwoDimension(&mySystem2, l0, N);
        myMeanValue.Matrix(mySystem2, l0, N, meanMat4);
        myfile4 << meanMat4 << endl;
    }
    wtime1 = omp_get_wtime() - wtime1;
    cout << "did MC 2+1 dims fixed step , Elapsed seconds = " << wtime1 << endl;

    //------------------------------------------------
    wtime2 = omp_get_wtime();
    for(int i = 0; i<TimeSteps; i++ ){
        myMonteCarlo.TwoDimensionGauss(&mySystem1, l0, N);     //these two for gauss
        myMeanValue.Matrix(mySystem1, l0, N, meanMat1);    //these two for gauss
        myfile1 << meanMat1 << endl;
    }
    wtime2 = omp_get_wtime() - wtime2;
    cout << "did MC 2+1 dims gauss step , Elapsed seconds = " << wtime2 << endl;

    //------------------------------------------------
    wtime3 = omp_get_wtime();
    for(int i = 0; i<TimeSteps; i++ ){
        myDiffSolver1.Explicit(meanMat2, 0.5*dt, l0, l0);  //diffsolver
        myfile2 << meanMat2 << endl;
    }
    wtime3 = omp_get_wtime() - wtime3;
    cout << "solved DiffEq 2+1 dims Explicit solver , Elapsed seconds = " << wtime3 << endl;

    //------------------------------------------------

    mat meanMat3 = zeros<mat>(int((1/l0))+1, int((1/l0))+1);
    wtime4 = omp_get_wtime();
    myDiffSolver2.Analytical(meanMat3, dt, l0, l0);  //diffsolver
    for(int i = 0; i<TimeSteps; i++ ){

        myfile3 << meanMat3 << endl;
    }
    wtime4 = omp_get_wtime() - wtime4;
    cout << "analytical Laplace, Elapsed seconds =  " << wtime4 << endl;

    //------------------------------------------------


    ofstream myfile8;
    myfile8.open("Implicit_jacobi.dat");
    double wtime8;
    wtime8 = omp_get_wtime();
    for(int i = 0; i<TimeSteps; i++ ){
        myDiffSolver3.ImplicitJacobi(meanMat5, int(1/l0), l0, dt);
        myfile8 << meanMat5 << endl;

    }
    wtime8 = omp_get_wtime() - wtime8;
    cout << "Implicit jacobi, Elapsed seconds =  " << wtime8 << endl;


    //close files
    myfile1.close(); myfile2.close(); myfile3.close();  myfile4.close(); myfile8.close();



    /* -------------------------------------------------
     * Solving Diffusion equation in 1+1 dimensions:
     *
     --------------------------------------------------*/




    mat u = zeros<mat>(length, TimeSteps); //use for non gauss
    vec meanVec = zeros(length);


    ofstream myfile5, myfile6;


    MonteCarloStep myMonteCarloStep3;
    MeanValue myMeanValue3;
    System mySystem3, mySystem4;

    double wtime5, wtime6;

      //adding particles to system (setting up initial state)
    for(int i = 0; i<N; i++){
        mySystem4.addParticle( 0.5 );
        mySystem3.addParticle( 0.5 );
    }




    //1d MC fixed step:
    wtime5 = omp_get_wtime();
       for(int i = 0; i<TimeSteps; i++){
        myMonteCarloStep3.OneDimension(&mySystem3, l0, N);
        myMeanValue3.Vector(mySystem3, l0, meanVec, N);
        u.col(i) = meanVec;
    }
       // Write to file
    myfile5.open("1d_MC_fixed.dat");
    for(int i = 0; i<length; i++){
        myfile5 << u.row(i);
        //cout << "there " << u.row(i) << endl;
    }
    myfile5 << endl;
    myfile5.close();
    wtime5 = omp_get_wtime() - wtime5;
    cout << "one dim MC fixed step, Elapsed seconds =  " << wtime5 << endl;

    //1d MC gauss step:
    wtime6 = omp_get_wtime();
       for(int i = 0; i<TimeSteps; i++){
        myMonteCarloStep3.OneDimensionGauss(&mySystem4, l0, N);
        myMeanValue3.Vector(mySystem4, l0, meanVec, N);
        u.col(i) = meanVec;
    }
       // Write to file
    myfile6.open("1d_MC_gauss.dat");
    for(int i = 0; i<length; i++){
        myfile6 << u.row(i);
        //cout << "there " << u.row(i) << endl;
    }
    myfile6 << endl;
    myfile6.close();
    wtime6 = omp_get_wtime() - wtime6;
    cout << "one dim MC gauss step, Elapsed seconds =  " << wtime6 << endl;


    cout << "finish!" << endl;
    return 0;
}

