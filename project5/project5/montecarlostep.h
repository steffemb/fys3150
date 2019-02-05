#ifndef MONTECARLOSTEP_H
#define MONTECARLOSTEP_H
#include <armadillo>
#include "system.h"



class MonteCarloStep
{
public:
    MonteCarloStep();

    void OneDimension(System *mySystem, double l0, int N);
    void TwoDimension(System *mySystem, double l0, int N);
    void RandNumTest();
    void TwoDimensionGauss(System *mySystem, double l0, int N);
    //void OneDimension(vec &xvector, double dt, int N, double l0, double D);

    void OneDimensionSimple(arma::vec &xvector, double dt, int N, double l0, double D, double TimeSteps, arma::vec &yvec, arma::mat &u);
    void OneDimensionGauss(System *mySystem, double l0, int N);
    void TwoDimension(arma::vec &xvector, arma::vec &yvector, double dt, int N, double l0, double D);
    //void MeanValue(int N, double l0, arma::vec &yvec, arma::vec xvector, arma::mat &u, double TimeSteps);
};
#endif // MONTECARLOSTEP_H
