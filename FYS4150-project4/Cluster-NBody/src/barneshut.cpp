#include<time.h>
#include "barneshut.h"

#include "octree.h"
#include "clustersystem.h"
#include "particle.h"

BarnesHut::BarnesHut(ClusterSystem *cs, Octree *oct){
    /*
     * Constructor for BarnesHut class.
     *
     * @param ClusterSystem, pointer to a system
     * @param Octree, pointer to octree
     */
    this->cs = cs;
    this->oct = oct;
    if ((this->cs == NULL) || (this->oct == NULL)){
        std::cout << "ClusterSystem or Octree NULL" << std::endl;
    }
    this->firstIndex = 0;
    this->lastIndex = 0;
}

void BarnesHut::calculateForces(){
    /*
     * Calculate forces on all particles in the octree.
     */
    this->cs->resetAllForces();

    int i;
#pragma omp parallel for num_threads(4)
    for (i = 0; i < this->cs->getNumberOfParticles(); i++){
        (this->oct)->calculateForceOnParticle( (this->cs)->getParticleAtIndex(i), this->tau , (this->cs)->getGravitationalConstant(), (this->cs)->getGravitationalSmoothing());
    }

    return;
}

ClusterSystem* BarnesHut::getClusterSystem(){
    return this->cs;
}

Octree* BarnesHut::getOctree(){
    return this->oct;
}

int BarnesHut::getIndexOfFirstParticle(){
    return this->firstIndex;
}

int BarnesHut::getIndexOfLastParticle(){
    return this->lastIndex;
}

void BarnesHut::setTau(double tau){
    /*
     * Set the BarnesHut parameter.
     *
     * @param double
     */
    this->tau = tau;
    return;
}

void BarnesHut::setClusterSystem(ClusterSystem *cs){
    this->cs = cs;
    return;
}

void BarnesHut::setOctree(Octree *oct){
    this->oct = oct;
    return;
}

void BarnesHut::setIndexOfFirstParticle(int index){
    this->firstIndex = index;
    return;
}

void BarnesHut::setIndexOfLastParticle(int index){
    this->lastIndex = index;
    return;
}
