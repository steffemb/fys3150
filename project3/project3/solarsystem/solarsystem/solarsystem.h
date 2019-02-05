#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include <celestialbody.h>
#include <vector>
#include <iostream>
using std::vector;
using namespace std;

class Solarsystem
{
public:
    vector<CelestialBody> objects;
    Solarsystem();
    void addCelestialBody(CelestialBody newObject);
    void Forces();
    void ForcesOldstyle(CelestialBody thisObject);
    double TotalEnergy();
    int NumberOfBodies();
    void resetForces();
    int lengthSquared(vec thisvector);
};

#endif // SOLARSYSTEM_H
