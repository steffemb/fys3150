#include "solarsystem.h"
#include <cmath>
#include <celestialbody.h>
#include <armadillo>


Solarsystem::Solarsystem()
{  
}

void Solarsystem::addCelestialBody(CelestialBody newObject) {
    objects.push_back(newObject);
    cout << "added object" << endl;
}

void Solarsystem::Forces() {

    resetForces();
    for(int i=0; i<NumberOfBodies(); i++) {
        //CelestialBody &body1 = objects[i];
        for(int j=i+1; j<NumberOfBodies(); j++) {
            //CelestialBody &body2 = objects[j];
            vec3 deltaRVector = objects[i].position - objects[j].position;
            double dr = sqrt(deltaRVector[0]*deltaRVector[0] + deltaRVector[1]*deltaRVector[1] + deltaRVector[2]*deltaRVector[2]);
            double G = 39.47; //8246681e2; //e17; //AU Msun^-1 (km/y)^2

            // Calculate the force and potential energy here
            double f = -(G * objects[i].mass * objects[j].mass) / (dr * dr * dr);

            objects[i].force += f*deltaRVector; // x, y, z -component of the force.
            objects[j].force -= f*deltaRVector; // Newtons third law.

            objects[i].potentialenergy += objects[i].mass*G*dr;

        }

        objects[i].kineticenergy += 0.5*objects[i].mass*lengthSquared(objects[i].velocity);
        //cout << "calculated forces " << objects[i].force << endl;
    }

}

int Solarsystem::NumberOfBodies()
{
    return objects.size();
}

double Solarsystem::TotalEnergy()
{
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    for(int i=0; i<NumberOfBodies(); i++) {
        potentialEnergy += objects[i].potentialenergy;
        kineticEnergy += objects[i].kineticenergy;
    }
    return kineticEnergy + potentialEnergy;
}

void Solarsystem::resetForces()
{
    for(int i=0; i < NumberOfBodies(); i++){
        objects[i].force = zeros(3);
        //cout << "reset forces = " << objects[i].force << endl;
    }
}

int Solarsystem::lengthSquared(vec thisvector)
{
    double total = 0;
    for(int i = 0; i<thisvector.size(); i++){
        total += thisvector[i]*thisvector[i];
    }
    return total;
}

