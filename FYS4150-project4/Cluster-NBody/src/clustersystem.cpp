#include "clustersystem.h"

#include<time.h>
#include<omp.h>
#include<armadillo>

#include "particle.h"
#include "barneshut.h"
#include "octree.h"
#include "gaussiandeviate.cpp"

using namespace arma;

// seeds for the random functions
//long IDUM_GAUSS = -412347645;
//long IDUM_UNIFORM = -134928459;
long IDUM_GAUSS = time(0);
long IDUM_UNIFORM = time(0) + 10;

ClusterSystem::ClusterSystem(){
    /*
     * Constructor for the ClusterSystem class.
     */
    this->numberOfParticles = 0;
    this->massParticle = 10.0;
}

ClusterSystem::ClusterSystem(int N, double systemSideLength){
    /*
     * Constructor for the ClusterSystem class.
     *
     * @param int, number of particles
     * @param double, size of system (diameter)
     */
    this->systemSideLength = systemSideLength;
    this->systemSideLengthOriginal = systemSideLength;
    this->massParticle = 10.0;
    this->numberOfParticles = 0;
    this->addANumberOfParticles(N);
    this->calculateGravitationalConstant();
}

void ClusterSystem::simulateBarnesHut(double t, double dt, bool useSmoothing, std::string filenameBase){
    /*
     * Simulate the system using BarnesHut method and Velocity Verlet.
     *
     * @param double simulation time
     * @param double step size
     * @param bool use smoothing
     * @param string base name for outputfiles
     */
    double time0, time1;
    time0 = omp_get_wtime();

    // calculate the smoothing factor
    this->gravitationSmoothing = 0;
    if (useSmoothing){
        this->gravitationSmoothing = 1.1 / pow(double(this->numberOfParticles), 0.28);
    }
    std::cout << "Smoothing factor: " << this->gravitationSmoothing << std::endl;

    this->timeCurrent = 0;
    // create octree and barneshut
    Octree oct(this, this->systemSideLength);
    BarnesHut bh(this, &oct);
    bh.setTau(0.5);
    bh.calculateForces();

    this->calculateEnergyBarnesHut();
    this->openOutputFiles(filenameBase);
    this->writeToOutputFiles(oct.getCenterOfMass());

    vec rNew = zeros<vec>(3);
    vec vNew = zeros<vec>(3);
    vec a = zeros<vec>(3);
    int iterations = 0;
    while (this->timeCurrent < t){
        this->timeCurrent += dt;

        // print progress
        if (int(100*this->timeCurrent / t) % 1 == 0){
            std::cout << "Progress: " << int(100 * this->timeCurrent / t) << " %\r";
            std::cout.flush();
        }

        // first part of verlet
        int i;
//#pragma omp parallel for
        for (i = 0; i < this->numberOfParticles; i++){
            Particle *currentParticle = this->particles[i];
            currentParticle->setEnergyPotential(0);
            a = currentParticle->getForce() / currentParticle->getMass();
            rNew = currentParticle->getPosition() + currentParticle->getVelocity() * dt + 0.5 * a * dt * dt;
            vNew = currentParticle->getVelocity() + 0.5 * a * dt;
            currentParticle->setNewPosition(rNew);
            currentParticle->setNewVelocity(vNew);

            // check wether the system has to be expanded due to bound particles outside the boundary.
            if (!currentParticle->getIsEjected()){
                if(fabs(max(rNew)) > 0.5*this->systemSideLength){
                    this->systemSideLength = 2*ceil(fabs(max(rNew)));
                }
                else if(fabs(min(rNew)) > 0.5*this->systemSideLength){
                    this->systemSideLength = 2*ceil(fabs(min(rNew)));
                }
            }
        }

        // recalculate forces
        oct.clear();
        oct.repopulateTree(this);
        bh.calculateForces();

        // last part of verlet
//#pragma omp parallel for
        for (i = 0; i < this->numberOfParticles; i++){
            Particle *currentParticle = this->particles[i];
            a = currentParticle->getForce() / currentParticle->getMass();
            vNew = currentParticle->getVelocity() + 0.5 * a * dt;
            currentParticle->setNewVelocity(vNew);
        }

        this->calculateEnergyBarnesHut();

        // do not write to file every iteration
        if (iterations == 24){
            iterations = 0;
            this->writeToOutputFiles(oct.getCenterOfMass());
        }
        iterations++;
    }
    this->closeOutputFiles();
    this->writeRadialDistances(oct.getCenterOfMass(), filenameBase);

    time1 = omp_get_wtime();
    std::cout << "Calculation time BarnesHut: " << (time1-time0) << " s." << std::endl;
    std::cout << "Time per iteration: " << (time1-time0)*(dt/t) << " s." << std::endl;

    int numberOfEjected = 0;
    int numberOfBound = 0;
    // count number of bound and ejected particles
    for (int i = 0; i < this->numberOfParticles; i++){
        if(this->particles[i]->getIsEjected()){
            numberOfEjected++;
        }
        if(this->particles[i]->getIsBound()){
            numberOfBound++;
        }
    }
    std::cout << "Ejected: " << numberOfEjected << " of " << this->numberOfParticles << std::endl;
    std::cout << "Bound: " << numberOfBound << " of " << this->numberOfParticles << std::endl;
    std::cout << "Final system size: " << this->systemSideLength << std::endl << std::endl;

    return;
}

void ClusterSystem::simulateNtoN(double t, double dt, bool useSmoothing, std::string filenameBase){
    /*
     * Simulate the system using direct method and Velocity Verlet.
     *
     * @param double simulation time
     * @param double step size
     * @param bool use smoothing
     * @param string base name for outputfiles
     */
    double time0, time1;
    time0 = omp_get_wtime();

    this->timeCurrent = 0;
    this->openOutputFiles(filenameBase);
    this->calculateEnergyNtoN();
    this->writeToOutputFiles();

    // calcualte the smoothing factor
    this->gravitationSmoothing = 0;
    if (useSmoothing){
        this->gravitationSmoothing = 1.1 / pow(double(this->numberOfParticles), 0.28);
    }
    std::cout << "Smoothing factor: " << this->gravitationSmoothing << std::endl;
    // find initial forces
    this->calculateForcesNtoN();

    vec rNew = zeros<vec>(3);
    vec vNew = zeros<vec>(3);
    vec a = zeros<vec>(3);
    int iterations = 0;
    while (this->timeCurrent < t){
        this->timeCurrent += dt;

        // print progress
        if (int(100*this->timeCurrent / t) % 1 == 0){
            std::cout << "Progress: " << int(100 * this->timeCurrent / t) << " %\r";
            std::cout.flush();
        }

        for (int i = 0; i < this->numberOfParticles; i++){
            Particle *currentParticle = this->particles[i];
            a = currentParticle->getForce() / currentParticle->getMass();
            rNew = currentParticle->getPosition() + currentParticle->getVelocity() * dt + 0.5 * a * dt * dt;
            vNew = currentParticle->getVelocity() + 0.5 * a * dt;
            currentParticle->setNewPosition(rNew);
            currentParticle->setNewVelocity(vNew);
        }

        this->calculateForcesNtoN();

        // last part of verlet
        for (int i = 0; i < this->numberOfParticles; i++){
            Particle *currentParticle = this->particles[i];
            a = currentParticle->getForce() / currentParticle->getMass();
            vNew = currentParticle->getVelocity() + 0.5 * a * dt;
            currentParticle->setNewVelocity(vNew);
        }

        this->calculateEnergyNtoN();

        // do not write data every step
        if (iterations == 24){
            iterations = 0;
            this->writeToOutputFiles();
        }
        iterations++;
    }

    this->closeOutputFiles();

    time1 = omp_get_wtime();
    std::cout << "Calculation time NtoN: " << (time1-time0) << " s." << std::endl;
    std::cout << "Time per iteration: " << (time1-time0)*(dt/t) << " s." << std::endl;

    int numberOfEjected = 0;
    int numberOfBound = 0;
    // count number of bound and ejected particles
    for (int i = 0; i < this->numberOfParticles; i++){
        if(this->particles[i]->getIsEjected()){
            numberOfEjected++;
        }
        if(this->particles[i]->getIsBound()){
            numberOfBound++;
        }
    }
    std::cout << "Ejected: " << numberOfEjected << " of " << this->numberOfParticles << std::endl;
    std::cout << "Bound: " << numberOfBound << " of " << this->numberOfParticles << std::endl;

    return;
}

void ClusterSystem::simulateNtoNRK4(double t, double dt, bool useSmoothing, std::string filenameBase){
    /*
     * Simulate the system using direct method and Runge-Kutta 4.
     *
     * @param double simulation time
     * @param double step size
     * @param bool use smoothing
     * @param string base name for outputfiles
     */
    double time0, time1;
    time0 = omp_get_wtime();

    this->timeCurrent = 0;
    this->openOutputFiles(filenameBase);
    this->calculateEnergyNtoN();
    this->writeToOutputFiles();

    // calcualte the smoothing factor
    this->gravitationSmoothing = 0;
    if (useSmoothing){
        this->gravitationSmoothing = 1.1 / pow(double(this->numberOfParticles), 0.28);
    }
    std::cout << "Smoothing factor: " << this->gravitationSmoothing << std::endl;


    int iterations = 0;
    while (this->timeCurrent < t){
        this->timeCurrent += dt;

        // print progress
        if (int(100*this->timeCurrent / t) % 1 == 0){
            std::cout << "Progress: " << int(100 * this->timeCurrent / t) << " %\r";
            std::cout.flush();
        }

        this->calculateForcesNtoN();
        this->RK4iterate(dt);

        this->calculateEnergyNtoN();

        // do not write to file every step
        if (iterations == 24){
            iterations = 0;
            this->writeToOutputFiles();
        }
        iterations++;
    }

    this->closeOutputFiles();

    time1 = omp_get_wtime();
    std::cout << "Calculation time NtoN: " << (time1-time0) << " s." << std::endl;
    std::cout << "Time per iteration: " << (time1-time0)*(dt/t) << " s." << std::endl;

    int numberOfEjected = 0;
    int numberOfBound = 0;
    // count number of bound and ejected particles
    for (int i = 0; i < this->numberOfParticles; i++){
        if(this->particles[i]->getIsEjected()){
            numberOfEjected++;
        }
        if(this->particles[i]->getIsBound()){
            numberOfBound++;
        }
    }
    std::cout << "Ejected: " << numberOfEjected << " of " << this->numberOfParticles << std::endl;
    std::cout << "Bound: " << numberOfBound << " of " << this->numberOfParticles << std::endl;

    return;
}

void ClusterSystem::RK4iterate(double dt){
    /*
     * Advance one step with Runge-Kutta 4.
     *
     * @param double dt, step size
     */
    int n = this->numberOfParticles;

    // matrices for storing intermediate step values
    mat k1 = zeros<mat>(3, n);
    mat k2 = zeros<mat>(3, n);
    mat k3 = zeros<mat>(3, n);
    mat k4 = zeros<mat>(3, n);
    mat l1 = zeros<mat>(3, n);
    mat l2 = zeros<mat>(3, n);
    mat l3 = zeros<mat>(3, n);
    mat l4 = zeros<mat>(3, n);

    // storing velocit and position at beginning of step
    mat currentVelocity = zeros<mat>(3, n);
    mat currentPosition = zeros<mat>(3, n);

    // calculate l1 and k1
    for (int i = 0; i < n; i++){
        Particle *currentParticle = this->particles[i];
        currentPosition.col(i) = currentParticle->getPosition();
        currentVelocity.col(i) = currentParticle->getVelocity();

        l1.col(i) = dt * currentParticle->getForce() / currentParticle->getMass();
        k1.col(i) = dt * currentParticle->getVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l1.col(i) / 2.0;
        vec newPos = currentPosition.col(i) + k1.col(i) / 2.0;
        currentParticle->setNewVelocity(newVel);
        currentParticle->setNewPosition(newPos);
    }

    this->calculateForcesNtoN();

    // calculate l2 and k2
    for (int i = 0; i < n; i++){
        Particle *currentParticle = this->particles[i];
        l2.col(i) = dt * currentParticle->getForce() / currentParticle->getMass();
        k2.col(i) = dt * currentParticle->getVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l2.col(i) / 2.0;
        vec newPos = currentPosition.col(i) + k2.col(i) / 2.0;
        currentParticle->setNewVelocity(newVel);
        currentParticle->setNewPosition(newPos);
    }

    this->calculateForcesNtoN();

    // calculate l3 and k3
    for (int i = 0; i < n; i++){
        Particle *currentParticle = this->particles[i];
        l3.col(i) = dt * currentParticle->getForce() / currentParticle->getMass();
        k3.col(i) = dt * currentParticle->getVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l3.col(i);
        vec newPos = currentPosition.col(i) + k3.col(i);
        currentParticle->setNewVelocity(newVel);
        currentParticle->setNewPosition(newPos);
    }

    this->calculateForcesNtoN();

    // calculate l4 and k4
    for (int i = 0; i < n; i++){
        Particle *currentParticle = this->particles[i];
        l4.col(i) = dt * currentParticle->getForce() / currentParticle->getMass();
        k4.col(i) = dt * currentParticle->getVelocity();

        // set final new velocity and position
        vec newVel = currentVelocity.col(i) + (l1.col(i) + 2*l2.col(i) + 2*l3.col(i) + l4.col(i)) / 6.0;
        vec newPos = currentPosition.col(i) + (k1.col(i) + 2*k2.col(i) + 2*k3.col(i) + k4.col(i)) / 6.0;
        currentParticle->setNewVelocity(newVel);
        currentParticle->setNewPosition(newPos);
    }
    return;
}

void ClusterSystem::addANumberOfParticles(int N){
    /*
     * Add particles to the system.
     * Uniform distributed in the spherical system.
     * Gaussian distributed mass.
     *
     * @param int N
     */
    this->numberOfParticles += N;

    double phi, r, theta, x, y, z, mass;
    for (int i = 0; i < N; i++){
        phi = 2 * M_PI * ran2(&IDUM_UNIFORM);
        r = 0.5 * this->systemSideLength * std::pow(ran2(&IDUM_UNIFORM), 1/3.0);
        theta = acos(1 - 2 * ran2(&IDUM_UNIFORM));
        x = r * sin(theta) * cos(phi);
        y = r * sin(theta) * sin(phi);
        z = r * cos(theta);
        mass = this->massParticle + gaussian_deviate(&IDUM_GAUSS);

        this->addParticleNoVelocity(x, y, z, mass);
    }
    return;
}

int ClusterSystem::getNumberOfParticles(){
    return this->numberOfParticles;
}

Particle* ClusterSystem::getParticleAtIndex(int index){
    if (index > this->numberOfParticles-1){
        return NULL;
    }
    return this->particles[index];
}

double ClusterSystem::getSystemSideLength(){
    return this->systemSideLength;
}

void ClusterSystem::calculateEnergyBarnesHut(){
    /*
     * Calculate the energy of the system when using BarnesHut.
     */
    this->energyKinetic = 0;
    this->energyPotential = 0;

    vec pos;
    for (int i = 0; i < this->numberOfParticles; i++){
        Particle *currentParticle = this->particles[i];
        // skip ejected particle
        if (currentParticle->getIsEjected()){
            continue;
        }

        pos = currentParticle->getPosition();
        currentParticle->calculateEnergyKinetic();

        // particle not bound
        if (currentParticle->getEnergyTotal() > 0){
            currentParticle->setIsBound(false);
            // particle unbound and outside system
            if (norm(pos) > 0.5*this->systemSideLengthOriginal){
                currentParticle->setIsEjected(true);
                continue;
            }
        }
        // insert bound test here, or else
        else{
            this->energyKinetic += currentParticle->getEnergyKinetic();
            this->energyPotential += 0.5* currentParticle->getEnergyPotential();
        }
    }
    return;
}

void ClusterSystem::calculateEnergyNtoN(){
    /*
     * Calculate the energy of the system.
     */
    this->energyKinetic = 0;
    this->energyPotential = 0;

    vec pos;
    double Ep;
    int i;
//#pragma omp parallel for private(Ep, pos) num_threads(4)
    for (i = 0; i < this->numberOfParticles; i++){
        Particle *currentParticle = this->particles[i];
        // skip ejected particles
        if (currentParticle->getIsEjected()){
            continue;
        }

        currentParticle->calculateEnergyKinetic();

        Ep = 0;
        for (int j = 0; j < this->numberOfParticles; j++){
            if (j == i){
                continue;
            }
            Particle *tempParticle = this->particles[j];
            if (tempParticle->getIsEjected()){
                continue;
            }
            pos = currentParticle->getPosition() - tempParticle->getPosition();
            Ep += - this->gravitationConstant * currentParticle->getMass() * tempParticle->getMass() / norm(pos);

        }
        currentParticle->setEnergyPotential(Ep);

        // check if particle is bound or ejected
        if (currentParticle->getEnergyTotal() > 0){
            currentParticle->setIsBound(false);
            // particle not bound and outside the system
            if (abs(norm(currentParticle->getPosition())) > this->systemSideLengthOriginal){
                currentParticle->setIsEjected(true);
                continue;
            }
            continue;
        }

        // increase the systems energy
        this->energyKinetic += currentParticle->getEnergyKinetic();
        this->energyPotential += 0.5 * currentParticle->getEnergyPotential();
    }
    return;
}

void ClusterSystem::calculateForcesNtoN(){
    /*
     * Directly calculate forces on all particles.
     */
    this->resetAllForces();
    double rPow;
    double smoothing2 = this->gravitationSmoothing * this->gravitationSmoothing;
    int i;
    vec r;
    double GMM;
    vec force;
#pragma omp parallel for private(r, GMM, rPow, force) //num_threads(4)
    for (i = 0; i < this->numberOfParticles; i++){
        Particle *currentParticle = this->particles[i];
        // skip ejected particles
        if (currentParticle->getIsEjected()){
            continue;
        }
        for (int j = i+1; j < this->numberOfParticles; j++){
            // skip ejected particles
            if (this->particles[j]->getIsEjected()){
                continue;
            }
            r = (currentParticle->getPosition() - this->particles[j]->getPosition());
            GMM = this->gravitationConstant * currentParticle->getMass() * this->particles[j]->getMass();
            rPow = norm(r);
            force = - normalise(r) * GMM / (rPow * rPow + smoothing2);
            currentParticle->increaseForce(force);
            this->particles[j]->increaseForce(-force);

        }
    }
    return;
}

void ClusterSystem::resetAllForces(){
    /*
     * Reset the forces for each particle.
     */
    for (int i = 0; i < this->numberOfParticles; i++){
        this->particles[i]->resetForce();
    }
    return;
}

void ClusterSystem::addParticle(vec position, vec velocity, double mass){
    Particle *newParticle = new Particle(position, velocity, mass);
    this->particles.push_back(newParticle);
    return;
}

void ClusterSystem::addParticleNoVelocity(double x, double y, double z, double mass){
    vec vel = zeros<vec>(3);
    vec pos = zeros<vec>(3);
    pos(0) = x; pos(1) = y; pos(2) = z;

    Particle *newParticle = new Particle(pos, vel, mass);
    this->particles.push_back(newParticle);

    return;
}

void ClusterSystem::calculateGravitationalConstant(){
    /*
     * Calculate the gravitational constant.
     */
    this->gravitationConstant = (M_PI * M_PI * pow(this->systemSideLength / 2.0, 3)) / (8 * this->numberOfParticles * this->massParticle);
    return;
}

double ClusterSystem::getGravitationalConstant(){
    return this->gravitationConstant;
}

double ClusterSystem::getGravitationalSmoothing(){
    return this->gravitationSmoothing;
}

void ClusterSystem::openOutputFiles(std::string filenameBase){
    std::string filenamePos = "positions" + filenameBase + ".xyz";
    std::string filenameEnergy = "energy" + filenameBase + ".dat";
    this->outputPosition = new std::ofstream(filenamePos.c_str());
    this->outputEnergy = new std::ofstream(filenameEnergy.c_str());
    return;
}

void ClusterSystem::writeToOutputFiles(){
    /*
     * Write data to files.
     */
    *outputPosition << this->numberOfParticles << std::endl;
    *outputPosition << "Optional line yo" << std::endl;
    for (int i = 0; i < this->numberOfParticles; i++){
        *outputPosition << "Ar " << this->particles[i]->getPositionX() << "\t" << this->particles[i]->getPositionY() << "\t" << this->particles[i]->getPositionZ() << "\t" << std::endl;
    }
    *outputEnergy << this->timeCurrent << "\t" << this->energyKinetic << "\t" << this->energyPotential << "\t" << std::endl;
    return;
}

void ClusterSystem::writeToOutputFiles(vec CoM){
    /*
     * Write data to files.
     *
     * @param vec center of mass of system
     */
    *outputPosition << this->numberOfParticles+1 << std::endl;
    *outputPosition << "Optional line yo" << std::endl;
    *outputPosition << "Fe " << CoM(0) << "\t" << CoM(1) << "\t" << CoM(2) << std::endl;
    for (int i = 0; i < this->numberOfParticles; i++){
        *outputPosition << "Ar " << this->particles[i]->getPositionX() << "\t" << this->particles[i]->getPositionY() << "\t" << this->particles[i]->getPositionZ() << "\t" << std::endl;
    }

    int numberOfBound = 0;
    int numberOfEjected = 0;
    for (int i = 0; i < this->numberOfParticles; i++){
        if (this->particles[i]->getIsBound()){
            numberOfBound++;
            continue;
        }
        if (this->particles[i]->getIsEjected()){
            numberOfEjected++;
        }
    }
    *outputEnergy << this->timeCurrent << "\t" << this->energyKinetic << "\t" << this->energyPotential << "\t" << numberOfBound << "\t" << numberOfEjected << std::endl;
    //*outputEnergy << this->timeCurrent << "\t" << this->energyKinetic << "\t" << this->energyPotential << "\t" << std::endl;
    return;
}

void ClusterSystem::closeOutputFiles(){
    this->outputPosition->close();
    this->outputEnergy->close();
    return;
}

void ClusterSystem::writeRadialDistances(vec CoM, std::string filenameBase){
    /*
     * Write positions for use in radial distribution.
     *
     * @param vec center of mass
     * @param string filenameBase
     */
    std::string filenameRad = "radial" + filenameBase + ".dat";
    std::ofstream outputRadial;
    outputRadial.open(filenameRad.c_str());

    outputRadial << CoM(0) << "\t" << CoM(1) << "\t" << CoM(2) << std::endl;

    for (int i = 0; i < this->numberOfParticles; i++){
        Particle *currentParticle = this->particles[i];
        if (currentParticle->getIsBound()){
            //outputRadial << norm(currentParticle->getPosition() - CoM) << std::endl;
            outputRadial << currentParticle->getPositionX() << "\t" << currentParticle->getPositionY() << "\t" << currentParticle->getPositionZ() << "\t" << std::endl;
        }
    }

    outputRadial.close();
    return;
}
