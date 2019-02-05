#include <iostream>
#include <solarsystem.h>
#include <armadillo>
#include <cmath>
#include <integrator.h>


using namespace arma;
using namespace std;

int main()
{
    double T = 10;
    double dt = 0.01;
    int steps = T/dt;

    Solarsystem mySolarsystem;


    double kk = 0.333;
    double x = kk*5.2;
    double y = sqrt(5.2*5.2 - x*x);
    //cout << sqrt(x*x + y*y) << endl;
    double vy = kk*3.25;
    double vx = sqrt(3.25*3.25 - vy*vy);
    double massJupiter = (1.9e-3);
    double CenterOfMassR = 1/(1 + 3e-6 + massJupiter)*(1.0*3e-6 + 5.2*massJupiter);
    double sunVeloc = sqrt(CenterOfMassR*0.0029);
    double jupiterVeloc = sqrt((5.2-CenterOfMassR)*0.0028/massJupiter);
    CelestialBody sun(-CenterOfMassR, 0, 0, 0, -sunVeloc, 0, 1.0);
    CelestialBody earth(1.0-CenterOfMassR, 0, 0, 0, 2*M_PI, 0, 3e-6);
    CelestialBody jupiter(5.2-CenterOfMassR, 0, 0, 0, jupiterVeloc, 0, massJupiter); //0.000954265748

    mySolarsystem.addCelestialBody(sun);
    mySolarsystem.addCelestialBody(earth);
    mySolarsystem.addCelestialBody(jupiter);
    mySolarsystem.Forces();
    //cout << jupiterVeloc << endl;
    Integrator myIntegrator;
    ofstream myfile;
    myfile.open("position_earth_euler.dat");
    for(int timestep = 0; timestep < steps + 1; timestep++) {
        // chose on one of the latter:
        myIntegrator.RK4(mySolarsystem, dt);
        //myIntegrator.Verlet(mySolarsystem, dt);
        //myIntegrator.EulerChromer(mySolarsystem, dt);

        cout << mySolarsystem.TotalEnergy() << endl; // to check energy conservation

        for(int i=0; i<mySolarsystem.NumberOfBodies(); i++) {
            myfile << float(mySolarsystem.objects[i].position[0]) << " " << float(mySolarsystem.objects[i].position[1]) << " ";
        }
        myfile << endl;
    }
    myfile.close();

    cout << "Earth position: " << mySolarsystem.objects[1].position << endl;
    return 0;
}

