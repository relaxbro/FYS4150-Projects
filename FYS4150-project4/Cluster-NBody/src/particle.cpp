#include "particle.h"

#include<armadillo>

using namespace arma;

Particle::Particle(vec position, vec velocity, double mass){
    /*
     * Constructor for particle class.
     *
     * @param vec position
     * @param vec velocity
     * @param double mass
     */
    this->position = position;
    this->velocity = velocity;
    this->mass = mass;
    this->force = zeros<vec>(3);
    this->energyKinetic = 0;
    this-> energyPotential = 0;
    this->isBound = true;
    this->isEjected = false;
}

vec Particle::getPosition(){
    return this->position;
}

double Particle::getPositionX(){
    return this->position(0);
}

double Particle::getPositionY(){
    return this->position(1);
}

double Particle::getPositionZ(){
    return this->position(2);
}

vec Particle::getVelocity(){
    return this->velocity;
}

double Particle::getMass(){
    return this->mass;
}

vec Particle::getForce(){
    return this->force;
}

double Particle::getForceX(){
    return this->force(0);
}

double Particle::getForceY(){
    return this->force(1);
}

double Particle::getForceZ(){
    return this->force(2);
}

double Particle::getEnergyTotal(){
    return this->energyKinetic + this->energyPotential;
}

double Particle::getEnergyKinetic(){
    return this->energyKinetic;
}

void Particle::calculateEnergyKinetic(){
    double v = norm(this->velocity);
    // calculate the kinetic energy
    this->energyKinetic = 0.5 * this->mass * v * v;
    return;
}

double Particle::getEnergyPotential(){
    return this->energyPotential;
}

void Particle::setEnergyPotential(double Ep){
    this->energyPotential = Ep;
    return;
}

void Particle::increaseEnergyPotential(double deltaEp){
    /*
     * Increase the potential energy of particle.
     *
     * @param double deltaEp
     */
    this->energyPotential += deltaEp;
    return;
}

void Particle::resetForce(){
    this->force(0) = 0;
    this->force(1) = 0;
    this->force(2) = 0;
    return;
}

void Particle::increaseForce(vec deltaForce){
    /*
     * Increase the force working on the particle.
     *
     * @param vec deltaforce
     */
    this->force += deltaForce;
    return;
}

void Particle::setNewPosition(vec newPosition){
    this->position = newPosition;
    return;
}

void Particle::setNewVelocity(vec newVelocity){
    this->velocity = newVelocity;
    return;
}

bool Particle::getIsBound(){
    return this->isBound;
}

void Particle::setIsBound(bool bound){
    this->isBound = bound;
    return;
}

bool Particle::getIsEjected(){
    return this->isEjected;
}

void Particle::setIsEjected(bool ejected){
    this->isEjected = ejected;
    return;
}
