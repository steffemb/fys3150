#include "celestialbody.h"

CelestialBody::CelestialBody(vec initialPosition, vec initialVelocity, double mass_) {
    position = initialPosition;
    velocity = initialVelocity;
    force = zeros(3);
    mass = mass_;

    //resetForce();
}

CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass_) {
    position = zeros(3);
    velocity = zeros(3);
    force = zeros(3);
    position[0] = x;
    position[1] = y;
    position[2] = z;

    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    mass = mass_;
    double kineticenergy = 0;
    double potentialenergy = 0;

    //resetForce();
}

void CelestialBody::resetForce() {
    force.zeros();
}
