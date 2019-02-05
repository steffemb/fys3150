#include "system.h"
#include "particle.h"
#include <vector>
#include <algorithm>

System::System()
{
}


void System::addParticle( double ypos ) {
    Particle newParticle((0.1)/2,ypos);
    particles.push_back(newParticle);
    //cout << "added particle bitches!" << endl;
}


