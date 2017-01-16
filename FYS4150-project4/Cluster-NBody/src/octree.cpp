#include "octree.h"

#include<armadillo>

#include "clustersystem.h"
#include "particle.h"

using namespace arma;

static const int NOT_AN_OCTANT = 8;

Octree::Octree(vec origin, double length, Particle *particle){
    /*
     * Constructor for Octree.
     *
     * @param vec origin, center of tree
     * @param double legth of each side
     * @param Particle pointer to particle to be placed in octree
     */
    this->origin = origin;
    this->halfOfLength = length / 2.0;
    this->isParent = false;
    this->centerOfMass = zeros<vec>(3);
    this->mass = 0;
    this->numberOfParticlesInside = 0;
    this->particle = NULL;
    this->children = NULL;
    this->addParticle(particle);
}

Octree::Octree(double x, double y, double z, double length, Particle *particle){
    /*
     * Constructor for Octree.
     *
     * @param double x-coordinate of origin
     * @param double y-coordinate of origin
     * @param double z-coordinate of origin
     * @param double legth of each side
     * @param Particle pointer to particle to be placed in octree
     */
    vec position = zeros<vec>(3);
    position(0) = x;
    position(1) = y;
    position(2) = z;
    this->origin = position;
    this->halfOfLength = length / 2.0;
    this->isParent = false;
    this->centerOfMass = zeros<vec>(3);
    this->mass = 0;
    this->numberOfParticlesInside = 0;
    this->particle = NULL;
    this->children = NULL;
    this->addParticle(particle);
}

Octree::Octree(ClusterSystem *system, double length){
    /*
     * Constructor for Octree.
     *
     * @param ClusterSystem pointer to system
     * @param double length of each side
     */
    this->origin = zeros<vec>(3);
    this->halfOfLength = length / 2.0;
    this->isParent = false;
    this->particle = NULL;
    this->children = NULL;
    this->centerOfMass = zeros<vec>(3);
    this->mass = 0;
    this->numberOfParticlesInside = 0;
    // add each particle in the system
    for (int i = 0; i < system->getNumberOfParticles(); i++){
        Particle *currentParticle = system->getParticleAtIndex(i);
        // dont add ejected particles
        if (currentParticle->getIsEjected()){
            continue;
        }
        this->addParticle(currentParticle);
    }
}

void Octree::calculateForceOnParticle(Particle *currentParticle, double tau, double gravitationalConstant, double smoothing){
    /*
     * Calculate the force on a particle->
     *
     * @param Particle pointer to particle
     * @param double tau BarnesHut paramter
     * @param double gravitational constant
     * @param double gravitational smoothing
     */
    // dont calculate force from itself
    if (currentParticle == this->particle){
        return;
    }

    // node is not parent. find force from particle
    if (!this->isParent){
        if (this->particle == NULL){
            return;
        }
        vec r = currentParticle->getPosition() - this->particle->getPosition();
        double GMM = gravitationalConstant * this->particle->getMass() * currentParticle->getMass();
        double rPow = norm(r);
        currentParticle->increaseEnergyPotential(-GMM / rPow);
        vec force = - normalise(r) * GMM / (rPow * rPow + smoothing * smoothing);
        currentParticle->increaseForce(force);
        return;
    }

    // distance from particle to center of mass of tree
    vec r = currentParticle->getPosition() - this->centerOfMass;

    // node is a parent node, check wether we can use the BarnesHut approximation
    // if not go to children of tree
    if ((2 * this->halfOfLength / norm(r)) > tau){
        this->children[0]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[1]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[2]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[3]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[4]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[5]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[6]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        this->children[7]->calculateForceOnParticle(currentParticle, tau, gravitationalConstant, smoothing);
        return;
    }
    // using approximation
    else {
        double GMM = gravitationalConstant * this->mass * currentParticle->getMass();
        double rPow = norm(r);
        currentParticle->increaseEnergyPotential(-GMM / rPow);
        vec force = - normalise(r) * GMM / (rPow * rPow + smoothing * smoothing);
        currentParticle->increaseForce(force);
        return;
    }

    return;
}

void Octree::addParticle(Particle *particle){
    /*
     * Add particle to tree.
     *
     * @param Particle pointer to particle
     */
    if (particle == NULL){
        return;
    }

    // check if particle is outside octant, should not happen
    if (this->getOctantForParticle(particle) == NOT_AN_OCTANT){
        std::cout << "Particle outside octant. Position:" << std::endl;
        std::cout << particle->getPosition() << std::endl;
        std::cout << "Energy: " << particle->getEnergyTotal() << std::endl;
        std::cout << "IsEjected: " << particle->getIsEjected() << std::endl;
        std::cout << "IsBound: " << particle->getIsBound() << std::endl;
        std::exit(1);
        return;
    }

    // particle inside node, increase count and update mass and center of mass
    this->numberOfParticlesInside += 1;
    //this->mass += particle->getMass();
    this->updateMassAndCenterOfMass(particle);

    // if node empty and not a parent set particle
    if ((this->particle == NULL) && (!this->isParent)){
        this->particle = particle;
        return;
    }

    // not parent and node full. create children and add to correct octant
    if (!this->isParent){
        this->createChildren();

        // add particle to correct octant
        this->children[this->getOctantForParticle(particle)]->addParticle(particle);
        // add allready containing particle to correct octant
        this->children[this->getOctantForParticle(this->particle)]->addParticle(this->particle);

        this->particle = NULL;
        return;
    }

    // allready parent, add to correct octant
    this->children[this->getOctantForParticle(particle)]->addParticle(particle);
    return;
}

void Octree::clear(){
    /*
     * Clear the octree, inlcudes deleting it to free memory.
     */
    if (!this->isParent){
        return;
    }

    for (int i = 0; i < 8; i++){
        if (this->children[i] == NULL){
            continue;
        }
        this->children[i]->clear();
        delete this->children[i];
        this->children[i] = NULL;
    }
    this->isParent = false;
    return;
}

void Octree::clear(Octree *oct){
    //this->particle == NULL;
    if (!this->isParent){
        std::cout << "tester delete" << std::endl;
        delete oct;
        return;
    }

    for (int i = 0; i < 8; i++){
        if (this->children[i] == NULL){
            continue;
        }
        this->clear(this->children[i]);
        //this->children[i]->clear(true);
        //delete this->children[i];
        //this->children[i] = NULL;
    }

    delete oct;
    //this->isParent = false;
    return;
}

void Octree::repopulateTree(ClusterSystem *system){
    /*
     * Repopulate the tree with particles in the system
     *
     * @param ClusterSystem pointer to system
     */
    this->halfOfLength = system->getSystemSideLength() / 2.0;

    // reset mass and number of children
    this->mass = 0;
    this->numberOfParticlesInside = 0;
    this->children = NULL;
    // add all particles
    for (int i = 0; i < system->getNumberOfParticles(); i++){
        Particle *currentParticle = system->getParticleAtIndex(i);
        // skip particles that have been ejected
        if (currentParticle->getIsEjected()){
            continue;
        }
        this->addParticle(currentParticle);
    }
    return;
}

Octree* Octree::getChild(int index){
    return this->children[index];
}

Particle* Octree::getParticle(){
    return this->particle;
}

bool Octree::getIsParant(){
    return this->isParent;
}

vec Octree::getOrigin(){
    return this->origin;
}

vec Octree::getCenterOfMass(){
    return this->centerOfMass;
}

int Octree::getNumberOfParticlesInside(){
    return this->numberOfParticlesInside;
}

double Octree::getMass(){
    return this->mass;
}

void Octree::updateMassAndCenterOfMass(Particle *particle){
    /*
     * Update the mass and center of mass.
     *
     * @param Particle pointer to particle
     */
    double newMass = this->mass + particle->getMass();
    this->centerOfMass = (this->mass * this->centerOfMass +  particle->getMass() * particle->getPosition()) / newMass;
    this->mass = newMass;
    return;
}

void Octree::createChildren(){
    /*
     * Create children.
     */
    if (this->isParent){
        std::cout << "should not be here, deleting kids" << std::endl;
        for (int i = 0; i < 8; i++){
            delete this->children[i];
            this->children[i] = NULL;
        }
    }

    this->children = new Octree*[8];
    // find new coordinates for children
    double newHalfOfLength = this->halfOfLength / 2.0;
    double x1 = this->origin(0) - newHalfOfLength;
    double y1 = this->origin(1) - newHalfOfLength;
    double z1 = this->origin(2) - newHalfOfLength;
    double x2 = this->origin(0) + newHalfOfLength;
    double y2 = this->origin(1) + newHalfOfLength;
    double z2 = this->origin(2) + newHalfOfLength;

    this->children[0] = new Octree(x1, y1, z1, this->halfOfLength, NULL);
    this->children[1] = new Octree(x2, y1, z1, this->halfOfLength, NULL);
    this->children[2] = new Octree(x2, y2, z1, this->halfOfLength, NULL);
    this->children[3] = new Octree(x1, y2, z1, this->halfOfLength, NULL);
    this->children[4] = new Octree(x1, y1, z2, this->halfOfLength, NULL);
    this->children[5] = new Octree(x2, y1, z2, this->halfOfLength, NULL);
    this->children[6] = new Octree(x2, y2, z2, this->halfOfLength, NULL);
    this->children[7] = new Octree(x1, y2, z2, this->halfOfLength, NULL);

    this->isParent = true;

    return;
}

int Octree::getOctantForParticle(Particle *particle){
    /*
     * Find correct octant for particle.
     *
     * @param Particle pointer to particle
     */
    if (particle == NULL){
        return NOT_AN_OCTANT;
    }

    vec position = particle->getPosition();

    // particle is outside the tree
    if ((position(0) < (origin(0) - halfOfLength)) || (position(0) > (origin(0) + halfOfLength))){
        return NOT_AN_OCTANT;
    }
    else if ((position(1) < (origin(1) - halfOfLength)) || (position(1) > (origin(1) + halfOfLength))){  
        return NOT_AN_OCTANT;
    }
    else if ((position(2) < (origin(2) - halfOfLength)) || (position(2) > (origin(2) + halfOfLength))){
        return NOT_AN_OCTANT;
    }

    // check coordinates to find correct octant
    if (position(0) < origin(0)){
        if (position(1) < origin(1)){
            if (position(2) < origin(2)){
                return 0;
            }
            else{
                return 4;
            }
        }
        else{
            if (position(2) < origin(2)){
                return 3;
            }
            else{
                return 7;
            }
        }
    }
    else{
        if (position(1) < origin(1)){
            if (position(2) < origin(2)){
                return 1;
            }
            else{
                return 5;
            }
        }
        else{
            if (position(2) < origin(2)){
                return 2;
            }
            else{
                return 6;
            }
        }
    }

    // if all faile. should not happen.
    return NOT_AN_OCTANT;
}

bool Octree::isLeaf(){
    // check if node is a leaf node
    if (this->children[0] == NULL){
        return true;
    }
    return false;
}

