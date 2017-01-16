#ifndef OCTREE_H
#define OCTREE_H

#include<armadillo>


class ClusterSystem;
class Particle;

class Octree{

public:
    Octree(arma::vec, double, Particle*);
    Octree(double, double, double, double, Particle*);
    Octree(ClusterSystem*, double);

    void calculateForceOnParticle(Particle*, double, double, double);
    void addParticle(Particle*);

    void clear();
    void clear(Octree*);
    void repopulateTree(ClusterSystem*);

    Octree* getChild(int);
    Particle* getParticle();
    bool getIsParant();
    arma::vec getOrigin();
    int getNumberOfParticlesInside();
    arma::vec getCenterOfMass();
    double getMass();
private:
    arma::vec origin;
    double halfOfLength;
    int numberOfParticlesInside;
    double mass;
    arma::vec centerOfMass;
    Octree **children;
    Particle *particle;
    bool isParent;

    void updateMassAndCenterOfMass(Particle*);
    bool isLeaf();
    void createChildren();
    int getOctantForParticle(Particle*);
};

#endif // OCTREE_H
