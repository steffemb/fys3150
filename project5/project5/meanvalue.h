#ifndef MEANVALUE_H
#define MEANVALUE_H
#include <armadillo>
#include "system.h"

class MeanValue
{
public:
    MeanValue();
    //void Vector(int N, double l0, arma::vec &yvec, arma::vec xvector);
    void Vector(System mySystem, double l0, arma::vec &meanVec, int N);
    void Matrix(System mySystem, double l0, int N, arma::mat &meanMat);
    void VectorSimple(int N, double l0, arma::vec &yvec, arma::vec xvector, arma::mat &u, double TimeSteps);
    //void MatrixGauss(System mySystem, double l0, arma::vec &meanVec, int N, arma::mat &meanMat);


};

#endif // MEANVALUE_H
