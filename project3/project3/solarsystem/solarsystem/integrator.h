#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <armadillo>
#include <solarsystem.h>

class Integrator
{
public:
    vec derivatives;
    Integrator();
    void EulerChromer(Solarsystem &mySolarsystem, double dt);
    void RK4(Solarsystem &mySolarsystem, double dt);
    void Verlet(Solarsystem &mySolarsystem, double dt);
};

#endif // INTEGRATOR_H
