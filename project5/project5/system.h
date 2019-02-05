#ifndef SYSTEM_H
#define SYSTEM_H

#include "particle.h"
#include <vector>
#include <iostream>

using std::vector;
using namespace std;


class System
{
public:
    System();
    vector<Particle> particles;
    void addParticle(double ypos);
    //void deleteParticle(Particle &thisParticle);
};

#endif // SYSTEM_H
