#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace arma;

class Solver
{
public:
    Solver();
    void ForwardEuler(mat &u, double i, double j, double dt, double dx);
    void BackwardEuler(mat &u, double i, double j, double dt, double dx);
    void CrankNichols(mat &u, double i, double j, double dt, double dx);
//    void ClosedForm(mat &u, double i, double j, double dt, double dx);


};

#endif // SOLVER_H
