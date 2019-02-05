#include "meanvalue.h"
#include <armadillo>
#include <iostream>
#include "system.h"

using namespace std;

MeanValue::MeanValue()
{
}


void MeanValue::Matrix(System mySystem, double l0, int N, arma::mat &meanMat){
    int numBoxesX = ((int((1/l0)) ));
    int numBoxesY = ((int((1/l0)) ));

    for(int k = 0; k< mySystem.particles.size(); k++){ //check all particles
        int boxXIndex = mySystem.particles[k].x_position*numBoxesX;
        int boxYIndex = mySystem.particles[k].y_position*numBoxesY;


        meanMat(boxYIndex,boxXIndex) += 1;
    }
    //cout << "at i,j =" << i << " " << j << "q = " << q << endl;

    double someNumber = N/(int((1/l0)));
    meanMat *= 1/someNumber;//double(60*(q)/((N)/(sigma))); // maybe(i, j) ?
}

void MeanValue::Vector(System mySystem, double l0, arma::vec &meanVec, int N){
    // mean value of N particles at position x
    int numBoxes = ((int((1/l0)) ));


        double q = 0; // reset counter
        for(int j = 0; j<mySystem.particles.size(); j++){
        int boxXIndex = mySystem.particles[j].x_position*numBoxes;

        meanVec(boxXIndex) += 1;

    }
    meanVec *= double(1.0/N);
}



void MeanValue::VectorSimple(int N, double l0, arma::vec &yvec, arma::vec xvector, arma::mat &u, double TimeSteps){
    // mean value of N particles at position x
    for(int k = 0; k<TimeSteps; k++){
        //arma::vec xvector = x.col(k);
        double xlength = xvector.n_elem;

        for(int i = 0; i< (int((1/l0))+1); i++){

            double q = 0; // reset counter
            for(int j = 0; j<xlength; j++){
                if ((i*l0-0.0000001 < xvector[j]) && (xvector[j]< i*l0+0.0000001)){
                    q += 1;
                }
            }
            yvec[i] += (q/N);
        }
        //cout << yvec << endl;
        u.col(k) = yvec;
    }
}

