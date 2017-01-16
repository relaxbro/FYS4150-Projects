#include "solarsystemobject.h"


SolarSystemObject::SolarSystemObject(string name, double mass, vec position, vec velocity){
    /*
     * Constructor for solar system object
     *
     * @param string name, name of object
     * @param double mass, mass of object
     * @param vec position, initial position
     * @param vec velocity, initial velocity
     */

    this->name = name;
    this->mass = mass;
    this->currentPosition = position;
    this->currentVelocity = velocity;
    this->forceG = zeros<vec>(3);

    this->hasLockedPosition = false;
}

string SolarSystemObject::getName(){
    return this->name;
}

double SolarSystemObject::getMass(){
    return this->mass;
}

vec SolarSystemObject::getCurrentPosition(){
    return this->currentPosition;
}

vec SolarSystemObject::getCurrentVelocity(){
    return this->currentVelocity;
}

vec SolarSystemObject::getForceG(){
    return this->forceG;
}

double SolarSystemObject::getForceGx(){
    return this->forceG(0);
}

double SolarSystemObject::getForceGy(){
    return this->forceG(1);
}

double SolarSystemObject::getForceGz(){
    return this->forceG(2);
}

vec SolarSystemObject::getDistanceToObject(SolarSystemObject otherObject){
    /*
     * Return distance to other object.
     *
     * @param SolarSystemObject otherObject
     * @return vec distance to other object
     */
    return (this->currentPosition - otherObject.getCurrentPosition());
}

bool SolarSystemObject::getHasLockedPosition(){
    return this->hasLockedPosition;
}

void SolarSystemObject::setHasLockedPosition(bool locked){
    /*
     * Function to set locked position
     *
     * @param bool locked
     */
    this->hasLockedPosition = locked;
    return;
}

void SolarSystemObject::setNewPosition(vec newPosition){
    /*
     * Function to set new position if not locked.
     *
     * @param vec newPosition
     */
    if (!this->hasLockedPosition){
        this->currentPosition = newPosition;
    }
    return;
}

void SolarSystemObject::setNewPosition(double x, double y, double z){
    /*
     * Function to set new position if not locked.
     *
     * @param double x
     * @param double y
     * @param double z
     */
    if (!this->hasLockedPosition){
        this->currentPosition(0) = x;
        this->currentPosition(1) = y;
        this->currentPosition(2) = z;
    }
    return;
}

void SolarSystemObject::setNewVelocity(vec newVelocity){
    /*
     * Function to set new velocity if position not locked.
     *
     * @param vec newVelocity
     */
    if (!this->hasLockedPosition){
        this->currentVelocity = newVelocity;
    }
    return;
}

void SolarSystemObject::setNewVelocity(double dx, double dy, double dz){
    /*
     * Function to set new velocity if position not locked.
     *
     * @param double vx
     * @param double vy
     * @param double vz
     */
    if (!this->hasLockedPosition){
        this->currentVelocity(0) = dx;
        this->currentVelocity(1) = dy;
        this->currentVelocity(2) = dz;
    }
    return;
}

void SolarSystemObject::increaseForceG(vec force){
    /*
     * Function to increase the force.
     *
     * @param vec force
     */
    this->forceG += force;
    return;
}

void SolarSystemObject::resetForceG(){
    /*
     * Function to reset the forces.
     */
    this->forceG(0) = 0;
    this->forceG(1) = 0;
    this->forceG(2) = 0;
    return;
}

