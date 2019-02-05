#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H
#include <armadillo>
using namespace arma;

class CelestialBody
{
public:
    vec position;
    vec velocity;
    vec force;

    double mass;
    double kineticenergy;
    double potentialenergy;

    CelestialBody(vec initialPosition, vec initialVelocity, double mass);
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void resetForce();
};

#endif // CELESTIALBODY_H
