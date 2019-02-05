#include "integrator.h"
#include <solarsystem.h>

Integrator::Integrator()
{
}

void Integrator::EulerChromer(Solarsystem &mySolarsystem, double dt){
    vec derivatives = zeros(3);
    mySolarsystem.Forces();
    for(int i = 0; i < mySolarsystem.objects.size(); i++) {
        CelestialBody &thisBody = mySolarsystem.objects[i];
        derivatives = thisBody.force/thisBody.mass;

        for(int k = 0; k < 3; k++){
            thisBody.velocity[k] += derivatives[k]*dt;
            thisBody.position[k] += thisBody.velocity[k]*dt;
        }

    }
}

void Integrator::RK4(Solarsystem &mySolarsystem, double dt){

    Solarsystem k1 = mySolarsystem;
    Solarsystem k2 = mySolarsystem;
    Solarsystem k3 = mySolarsystem;
    Solarsystem k4 = mySolarsystem;

    k1.Forces(); // We have k1

    for(int i = 0; i < k2.objects.size(); i++) {
        CelestialBody &bodyInPreviousSystem = k1.objects[i];
        CelestialBody  &bodyInThisSystem = k2.objects[i];

        bodyInThisSystem.velocity += bodyInPreviousSystem.force/bodyInPreviousSystem.mass*dt/2;
        bodyInThisSystem.position += bodyInPreviousSystem.velocity*dt/2;
    }
    k2.Forces(); // We have k2

    for(int i = 0; i < k3.objects.size(); i++) {
        CelestialBody &bodyInPreviousSystem = k2.objects[i];
        CelestialBody &bodyInThisSystem = k3.objects[i];

        bodyInThisSystem.velocity += bodyInPreviousSystem.force/bodyInPreviousSystem.mass*dt/2;
        bodyInThisSystem.position += bodyInPreviousSystem.velocity*dt/2;
    }
    k3.Forces(); // We have k3

    for(int i = 0; i < k4.objects.size(); i++) {
        CelestialBody &bodyInPreviousSystem = k3.objects[i];
        CelestialBody &bodyInThisSystem = k4.objects[i];

        bodyInThisSystem.velocity += bodyInPreviousSystem.force/bodyInPreviousSystem.mass*dt;
        bodyInThisSystem.position += bodyInPreviousSystem.velocity*dt;
    }
    k4.Forces(); // We have k4


    for(int i = 0; i < mySolarsystem.objects.size(); i++) {
        CelestialBody &body = mySolarsystem.objects[i];
        CelestialBody &body1 = k1.objects[i];
        CelestialBody &body2 = k2.objects[i];
        CelestialBody &body3 = k3.objects[i];
        CelestialBody &body4 = k4.objects[i];

        vec velocity = 1.0/6*(body1.velocity + 2*body2.velocity + 2*body3.velocity + body4.velocity);
        vec acceleration = 1.0/6*(body1.force + 2*body2.force + 2*body3.force + body4.force)/body.mass;

        body.position += velocity*dt;
        body.velocity += acceleration*dt;
    }

}

void Integrator::Verlet(Solarsystem &mySolarsystem, double dt){
    //integrator witch applies Verlet leapfrog method
    vec derivatives = zeros(3);
    Solarsystem halfstepSolarsystem = mySolarsystem;
    mySolarsystem.Forces();
    for(int i = 0; i < mySolarsystem.objects.size(); i++) {
        CelestialBody &thisBody = mySolarsystem.objects[i];
        CelestialBody &halfstepBody = halfstepSolarsystem.objects[i];
        derivatives = thisBody.force/thisBody.mass;
        halfstepBody.velocity = thisBody.velocity + derivatives*(dt/2.0);
        thisBody.position = halfstepBody.position + dt*halfstepBody.velocity;
    }
    mySolarsystem.Forces();
    for(int i = 0; i < mySolarsystem.objects.size(); i++) {
        CelestialBody &thisBody = mySolarsystem.objects[i];
        CelestialBody &halfstepBody = halfstepSolarsystem.objects[i];
        derivatives = thisBody.force/thisBody.mass;
        thisBody.velocity = halfstepBody.velocity + (dt/2.0)*derivatives;

    }
}
