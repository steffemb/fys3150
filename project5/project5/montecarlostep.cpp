#include <random> //c++ 11
#include <iostream>
#include "montecarlostep.h"
#include <armadillo>
#include <cmath>
#include <omp.h>
#include "system.h"
#include <vector>
#include <algorithm>

using namespace arma;
using namespace std;

MonteCarloStep::MonteCarloStep()
{
}

void MonteCarloStep::RandNumTest(){
    std::random_device rd;                //c++ 11
    std::mt19937 gen(rd());               //c++ 11
    std::uniform_real_distribution<> dis(0, 1);         //c++ 11
    int n = 100000;
    double first, second, third, fourth;
    for(int i=0; i<n; i++){
        double RandomNumber = dis(gen);          //c++ 11
        cout << "randomnum: " << RandomNumber << endl;
        if(RandomNumber <=0.25 && RandomNumber >=0){
            first += 1;
        }
        if(RandomNumber <=0.5 && RandomNumber >0.25){
            second += 1;
        }
        if(RandomNumber <=0.75 && RandomNumber >0.5){
            third +=1;
        }
        if(RandomNumber <=1 && RandomNumber >0.75){
            fourth += 1;
        }

    }
    cout << "if these numbers are not equal u are shit out of luck: " << first << " " << second << " " << third << " " << fourth<< endl;
}

void MonteCarloStep::TwoDimensionGauss(System *mySystem, double l0, int N){

    //new style!

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1/sqrt(2));
    std::random_device rd;                //c++ 11
    std::mt19937 gen(rd());               //c++ 11
    std::uniform_real_distribution<> dis(0, 1);          //c++ 11
    int n = mySystem->particles.size(); //dynamic particle number
    double eps = 0.00001; //how close to the boundary conditions can particles be?
    double original_l0 = l0;
    int nAdd = 0;
    int numberOfParticlesInBox = N/int(1/l0);
    int numBoxesX = (1/original_l0);
    int numBoxesY = (1/original_l0);
    vec boxCounter = zeros(numBoxesY);
    for(int i = 1; i<numBoxesY-1; i++){
        boxCounter[i] = N/(int(1/l0));
    }
    //cout << boxCounter << endl;
    //cout << "total number of particles: " << n << endl;
    for(int i=0; i<n; i++){ //itteration over particles
        int boxIndexXBefore = mySystem->particles[i].x_position*numBoxesX;
        int boxIndexYBefore = mySystem->particles[i].y_position*numBoxesY;

        double RandomNumber = dis(gen);          //c++ 11
        double RandomNumberGauss = distribution(generator);       //c++ 11
        l0 = original_l0*RandomNumberGauss;


        //test------------
        if(RandomNumber <= 0.25){
            mySystem->particles[i].x_position += l0;
        }
        if(RandomNumber > 0.25 && RandomNumber <= 0.5){
            mySystem->particles[i].x_position += l0;
        }
        if(RandomNumber > 0.5 && RandomNumber <= 0.75){
            mySystem->particles[i].y_position += l0;
        }
        if(RandomNumber > 0.75 && RandomNumber <= 1.0){
            mySystem->particles[i].y_position += l0;
        }
        int boxIndexXAfter = mySystem->particles[i].x_position*numBoxesX;
        int boxIndexYAfter = mySystem->particles[i].y_position*numBoxesY;
        //cout <<"xpos: " <<mySystem->particles[i].x_position << " boxIndexXAfter: " << boxIndexXAfter  << endl;

        //cout << mySystem->particles[i].x_position << " " <<  boxIndexXBefore << " " <<  boxIndexXAfter  << " " << numBoxesX << " " << mySystem->particles[i].x_position*numBoxesX << endl;
        if(boxIndexXBefore == 0 && boxIndexXAfter != boxIndexXBefore && mySystem->particles[i].x_position > 0) {
            // Reduce counter in box (0, )
            boxCounter(boxIndexYBefore) -= 1;
        }

        if(boxIndexXAfter == 0 && boxIndexXBefore != boxIndexXAfter && mySystem->particles[i].x_position > 0) {
            // Reduce counter in box (0, )
            boxCounter(boxIndexYBefore) += 1;
        }


        bool removeParticle = false;

        //boundary conditions:
        if(mySystem->particles[i].y_position >= 1.0-original_l0) removeParticle = true;
        if(mySystem->particles[i].y_position <= 0.0+original_l0) removeParticle = true;
        if(mySystem->particles[i].x_position >= 1.0-original_l0) removeParticle = true;
        if(mySystem->particles[i].x_position <= 0) removeParticle = true;


        if(removeParticle == true){
            if(boxIndexXBefore == 0 ) {
                // Reduce counter in box (0, )
                boxCounter(boxIndexYBefore) -= 1;
            }
        }

        if(boxIndexXBefore == 0 && removeParticle == false && boxIndexXAfter == 0 && boxIndexYBefore != boxIndexYAfter){
            boxCounter(boxIndexYBefore) -= 1;
            boxCounter(boxIndexYAfter) += 1;
        }

        //cout << boxCounter << mySystem->particles[i].x_position << endl;
        if(removeParticle) {
            swap(mySystem->particles[i], mySystem->particles.back());
            mySystem->particles.pop_back();
            i--;
            n--;
            continue;
        }
        //        //hard reset at pos (0, 0.5)
        //        if((-eps < mySystem->particles[i].x_position) && (mySystem->particles[i].x_position< original_l0)){
        //            //original_l0 = original_l0/2.0;
        //            //if(( (original_l0/2.0)+0.5-eps < mySystem->particles[i].y_position) && (mySystem->particles[i].y_position< (original_l0/2.0)+0.5+eps)){
        //                mySystem->particles.erase(mySystem->particles.begin()+i);
        //                i--;
        //                n--;
        //                continue;
        //            //}
        //        }
        //        //cout << "xpos: " << mySystem->particles[i].x_position << "ypos: " << mySystem->particles[i].y_position << endl;
    }

    //    //add N particles to system
    //    for(int i = 0; i<N; i++){
    //        mySystem->addParticle(0.5);
    //    }

    //cout << boxCounter << endl;
    //count particles in first bin x
    //or add N/int(1/original_l0) particles for each bin
    for(int k = 1; k<(numBoxesX)-1; k++){
        int temp = k;
        //cout<< "must be positive: " << numberOfParticlesInBox-boxCounter(temp) << " boxCounter = " << boxCounter[temp] << endl;
        for(int i = 0; i<numberOfParticlesInBox-boxCounter((temp)); i++){
            mySystem->addParticle(original_l0*(temp)+(original_l0/2));
        }
    }
}


void MonteCarloStep::TwoDimension(System *mySystem, double l0, int N){

    //new style!

    std::random_device rd;                //c++ 11
    std::mt19937 gen(rd());               //c++ 11
    std::uniform_real_distribution<> dis(0, 1);          //c++ 11
    int n = mySystem->particles.size(); //dynamic particle number
    double original_l0 = l0;
    int numberOfParticlesInBox = N/int(1/l0);
    int numBoxesX = (1/original_l0);
    int numBoxesY = (1/original_l0);
    vec boxCounter = zeros(numBoxesY);
    for(int i = 1; i<numBoxesY-1; i++){
        boxCounter[i] = N/(int(1/l0));
    }
    //cout << boxCounter << endl;
    //cout << "total number of particles: " << n << endl;
    for(int i=0; i<n; i++){ //itteration over particles
        int boxIndexXBefore = mySystem->particles[i].x_position*numBoxesX;
        int boxIndexYBefore = mySystem->particles[i].y_position*numBoxesY;

        double RandomNumber = dis(gen);          //c++ 11

        //test------------
        if(RandomNumber <= 0.25){
            mySystem->particles[i].x_position += l0;
        }
        if(RandomNumber > 0.25 && RandomNumber <= 0.5){
            mySystem->particles[i].x_position -= l0;
        }
        if(RandomNumber > 0.5 && RandomNumber <= 0.75){
            mySystem->particles[i].y_position += l0;
        }
        if(RandomNumber > 0.75 && RandomNumber <= 1.0){
            mySystem->particles[i].y_position -= l0;
        }
        int boxIndexXAfter = mySystem->particles[i].x_position*numBoxesX;
        int boxIndexYAfter = mySystem->particles[i].y_position*numBoxesY;

        if(boxIndexXBefore == 0 && boxIndexXAfter != boxIndexXBefore && mySystem->particles[i].x_position > 0) {
            boxCounter(boxIndexYBefore) -= 1;
        }

        if(boxIndexXAfter == 0 && boxIndexXBefore != boxIndexXAfter && mySystem->particles[i].x_position > 0) {
            boxCounter(boxIndexYBefore) += 1;
        }


        bool removeParticle = false;

        //boundary conditions:
        if(mySystem->particles[i].y_position >= 1.0-original_l0) removeParticle = true;
        if(mySystem->particles[i].y_position <= 0.0+original_l0) removeParticle = true;
        if(mySystem->particles[i].x_position >= 1.0-original_l0) removeParticle = true;
        if(mySystem->particles[i].x_position <= 0) removeParticle = true;


        if(removeParticle == true){
            if(boxIndexXBefore == 0 ) {
                boxCounter(boxIndexYBefore) -= 1;
            }
        }

        if(boxIndexXBefore == 0 && removeParticle == false && boxIndexXAfter == 0 && boxIndexYBefore != boxIndexYAfter){
            boxCounter(boxIndexYBefore) -= 1;
            boxCounter(boxIndexYAfter) += 1;
        }

        //cout << boxCounter << mySystem->particles[i].x_position << endl;
        if(removeParticle) {
            swap(mySystem->particles[i], mySystem->particles.back());
            mySystem->particles.pop_back();
            i--;
            n--;
            continue;
        }
    }

    for(int k = 1; k<(numBoxesX)-1; k++){
        int temp = k;
        //cout<< "must be positive: " << numberOfParticlesInBox-boxCounter(temp) << " boxCounter = " << boxCounter[temp] << endl;
        for(int i = 0; i<numberOfParticlesInBox-boxCounter((temp)); i++){
            mySystem->addParticle(original_l0*(temp)+(original_l0/2));
        }
    }
}

void MonteCarloStep::OneDimension(System *mySystem, double l0, int N){



    //new style!

    std::random_device rd;                //c++ 11
    std::mt19937 gen(rd());               //c++ 11
    std::uniform_real_distribution<> dis(0, 1);          //c++ 11
    int n = mySystem->particles.size(); //dynamic particle number
    double original_l0 = l0;
    int numberOfParticlesInBox = N;
    int numBoxesX = (1/original_l0);
    int boxCounter = N;

    //cout << boxCounter << endl;
    //cout << "total number of particles: " << n << endl;
    for(int i=0; i<n; i++){ //itteration over particles
        int boxIndexXBefore = mySystem->particles[i].x_position*numBoxesX;


        double RandomNumber = dis(gen);          //c++ 11

        //test------------
        if(RandomNumber <= 0.5){
            mySystem->particles[i].x_position += l0;
        }
        if(RandomNumber > 0.5 && RandomNumber <= 1){
            mySystem->particles[i].x_position -= l0;
        }

        int boxIndexXAfter = mySystem->particles[i].x_position*numBoxesX;


        if(boxIndexXBefore == 0 && boxIndexXAfter != boxIndexXBefore && mySystem->particles[i].x_position > 0) {
            boxCounter -= 1;
        }

        if(boxIndexXAfter == 0 && boxIndexXBefore != boxIndexXAfter && mySystem->particles[i].x_position > 0) {
            boxCounter += 1;
        }


        bool removeParticle = false;

        //boundary conditions:
        if(mySystem->particles[i].x_position >= 1.0) removeParticle = true;
        if(mySystem->particles[i].x_position <= 0) removeParticle = true;


        if(removeParticle == true){
            if(boxIndexXBefore == 0 ) {
                boxCounter -= 1;
            }
        }

        //cout << boxCounter << mySystem->particles[i].x_position << endl;
        if(removeParticle) {
            swap(mySystem->particles[i], mySystem->particles.back());
            mySystem->particles.pop_back();
            i--;
            n--;
            continue;
        }
    }

    for(int i = 0; i<numberOfParticlesInBox-boxCounter; i++){
        mySystem->addParticle(0.5);
    }

}



void MonteCarloStep::OneDimensionGauss(System *mySystem, double l0, int N){

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1/sqrt(2));
    std::random_device rd;                                          //c++ 11
    std::mt19937 gen(rd());                                         //c++ 11
    std::uniform_real_distribution<> dis(0, 1);                     //c++ 11
    int n = mySystem->particles.size(); //dynamic particle number
    double original_l0 = l0;
    int numberOfParticlesInBox = N;
    int numBoxesX = (1/original_l0);
    int boxCounter = N;

    //cout << boxCounter << endl;
    //cout << "total number of particles: " << n << endl;


    for(int i=0; i<n; i++){ //itteration over particles
        int boxIndexXBefore = mySystem->particles[i].x_position*numBoxesX;


        double RandomNumber = dis(gen);                             //c++ 11
        double RandomNumberGauss = distribution(generator);         //c++ 11
        l0 = original_l0*RandomNumberGauss;                         //c++ 11


        if(RandomNumber <= 0.5){
            mySystem->particles[i].x_position += l0;
        }
        if(RandomNumber > 0.5 && RandomNumber <= 1){
            mySystem->particles[i].x_position -= l0;
        }

        int boxIndexXAfter = mySystem->particles[i].x_position*numBoxesX;


        if(boxIndexXBefore == 0 && boxIndexXAfter != boxIndexXBefore && mySystem->particles[i].x_position > 0) {
            boxCounter -= 1;
        }

        if(boxIndexXAfter == 0 && boxIndexXBefore != boxIndexXAfter && mySystem->particles[i].x_position > 0) {
            boxCounter += 1;
        }


        bool removeParticle = false;

        //boundary conditions:
        if(mySystem->particles[i].x_position >= 1.0) removeParticle = true;
        if(mySystem->particles[i].x_position <= 0) removeParticle = true;


        if(removeParticle == true){
            if(boxIndexXBefore == 0 ) {
                boxCounter -= 1;
            }
        }


        if(removeParticle) {
            swap(mySystem->particles[i], mySystem->particles.back());
            mySystem->particles.pop_back();
            i--;
            n--;
            continue;
        }
    }

    for(int i = 0; i<numberOfParticlesInBox-boxCounter; i++){
        mySystem->addParticle(0.5);
    }


}


void MonteCarloStep::OneDimensionSimple(vec &xvector, double dt, int N, double l0, double D, double TimeSteps, arma::vec &yvec, arma::mat &u){

    //This is only a test NOT finalized
    //
    //

    std::random_device rd;                //c++ 11
    std::mt19937 gen(rd());               //c++ 11
    std::uniform_real_distribution<> dis(0, 1);          //c++ 11
    //cout << xvector.n_elem << endl;
    for(int k = 0; k<TimeSteps; k++){

        int numWalkers = xvector.n_elem;
        for (int i=0; i<numWalkers; i++){
            double RandomNumber = dis(gen);          //c++ 11

            //double RandomNumber = (rand()/(double)(RAND_MAX));   //c++ 98
            //cout << "RandomNumber " << RandomNumber << endl;
            if ((-0.0000001 < xvector[i]) && (xvector[i]< 0.0000001)){
                //cout << "to iiS!" << endl;
                if (RandomNumber <= 0.5){
                    xvector[i] += l0;
                    //cout << "jepps!" << xvector.n_elem << endl;
                }
                else if (RandomNumber > 0.5){
                    xvector[i] = 0;
                }
                continue;
            }
            else {
                if (RandomNumber <= 0.5){
                    xvector[i] += l0;
                }
                else if (RandomNumber > 0.5){
                    xvector[i] -= l0;
                }

            }

            if (1.0-0.0000001 <= xvector[i]) {
                xvector[i] = 0;
                std::swap(xvector[i], xvector[numWalkers-1]);
                xvector.resize(xvector.n_elem -1);
                numWalkers--;
                i--;
                // xvector[i].erase();
            }
        }
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
        //myMeanvalue(N, l0, yvec, xvector, u, TimeSteps);
        u.col(k) = yvec;
    }
}




