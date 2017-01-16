#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <math.h>
#include <vector>
#include <armadillo>
#include <boost/algorithm/string.hpp>

#include "solarsystemobject.h"
#include "verlet.h"
#include "rungekutta4.h"

using namespace std;
using namespace arma;

class SolarSystem{

public:
    SolarSystem(int, string);

    int getDimension();
    int getNumberOfCelestialBodies();

    void createCelestialBody(string, bool = false);
    void simulateVerlet(double, double, string);
    void simulateRK4(double, double, string);

private:
    //SolarSystem(const SolarSystem&);

    int dimension;
    double simulationTime;
    double simulationTimeCurrent;
    double simulationTimeStep;
    int simulationSteps;
    string bodiesFilename;
    int numberOfIntegrationSteps;

    double gravitationalConstant;

    double massSystem;
    double currentEnergyPotential;
    double currentEnergyKinetic;
    vec currentAngulerMomentum;

    vector<SolarSystemObject> celestialBodies;

    ofstream* outfilePosition;
    ofstream* outfileVelocity;
    ofstream* outfileForce;
    ofstream* outfileEnergy;
    ofstream* outfileAngularMomentum;

    void addCelestialBody(SolarSystemObject);

    void findForceOnBody(int);
    void findForceOnAllBodies();
    void resetForceOnAllBodies();
    void findEnergyKinetic();
    void findEnergyPotential();
    void resetEnergy();
    void findAngularMomentum();
    void resetAngularMomentum();

    void RK4iterate(double);

    void openOutputFiles(string);
    void writeOutputFiles();
    void closeOutputFiles();

};

#endif // SOLARSYSTEM_H
