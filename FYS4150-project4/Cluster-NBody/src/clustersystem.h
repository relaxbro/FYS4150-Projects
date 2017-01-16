#ifndef CLUSTERSYSTEM_H
#define CLUSTERSYSTEM_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<armadillo>
#include<math.h>


class Particle;
class BarnesHut;
class Octree;

class ClusterSystem{

public:
    ClusterSystem();
    ClusterSystem(int, double);

    void addANumberOfParticles(int);

    int getNumberOfParticles();
    double getGravitationalConstant();
    double getGravitationalSmoothing();
    double getSystemSideLength();
    Particle* getParticleAtIndex(int);
    void resetAllForces();

    void simulateBarnesHut(double, double, bool, std::string);
    void simulateNtoN(double, double, bool, std::string);
    void simulateNtoNRK4(double, double , bool, std::string);
    void RK4iterate(double);

private:
    int numberOfParticles;
    double systemSideLength;
    double systemSideLengthOriginal;
    double massParticle;
    double gravitationConstant;
    double gravitationSmoothing;
    double timeCurrent;

    double energyKinetic;
    double energyPotential;

    std::vector<Particle*> particles;

    void addParticle(arma::vec, arma::vec, double);
    void addParticleNoVelocity(double, double, double, double);

    void calculateGravitationalConstant();
    void calculateEnergyBarnesHut();
    void calculateEnergyNtoN();
    void calculateForcesNtoN();

    void openOutputFiles(std::string);
    void writeToOutputFiles();
    void writeToOutputFiles(arma::vec);
    void closeOutputFiles();
    void writeRadialDistances(arma::vec, std::string);

    std::ofstream *outputPosition;
    std::ofstream *outputEnergy;
};

#endif // CLUSTERSYSTEM_H
