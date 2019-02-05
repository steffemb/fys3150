#ifndef NUMERICALDIFF_H
#define NUMERICALDIFF_H

#include <armadillo>

class NumericalDiff
{
public:
    NumericalDiff();
    void Explicit(arma::mat &u, double dt, double dx, double dy);
    void Analytical(arma::mat &u, double dt, double dx, double dy);
    void ImplicitJacobi(arma::mat &u, int length, double h, double dt);

};

#endif // NUMERICALDIFF_H
