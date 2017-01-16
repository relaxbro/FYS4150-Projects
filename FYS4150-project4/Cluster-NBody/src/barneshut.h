#ifndef BARNESHUT_H
#define BARNESHUT_H


class Octree;
class ClusterSystem;
class Particle;

class BarnesHut{

public:
    BarnesHut(ClusterSystem*, Octree*);

    void calculateForces();

    ClusterSystem* getClusterSystem();
    Octree* getOctree();
    int getIndexOfFirstParticle();
    int getIndexOfLastParticle();

    void setTau(double);
    void setClusterSystem(ClusterSystem*);
    void setOctree(Octree*);
    void setIndexOfFirstParticle(int);
    void setIndexOfLastParticle(int);
private:
    ClusterSystem *cs;
    Octree *oct;

    int firstIndex;
    int lastIndex;
    double tau;
};

#endif // BARNESHUT_H
