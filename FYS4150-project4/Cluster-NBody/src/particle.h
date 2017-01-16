#ifndef PARTICLE_H
#define PARTICLE_H

#include<armadillo>


class Particle{

public:
    Particle(arma::vec, arma::vec, double);

    //inline arma::vec getPosition() {return position;}
    arma::vec getPosition();
    double getPositionX();
    double getPositionY();
    double getPositionZ();
    arma::vec getVelocity();
    double getMass();
    arma::vec getForce();
    double getForceX();
    double getForceY();
    double getForceZ();
    double getEnergyTotal();
    double getEnergyKinetic();
    void calculateEnergyKinetic();
    double getEnergyPotential();
    void setEnergyPotential(double);
    void increaseEnergyPotential(double);
    void resetForce();
    void increaseForce(arma::vec);
    void setNewPosition(arma::vec);
    void setNewVelocity(arma::vec);

    bool getIsBound();
    void setIsBound(bool);
    bool getIsEjected();
    void setIsEjected(bool);
private:
    double mass;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;

    arma::vec position;
    arma::vec velocity;
    arma::vec force;
    double energyKinetic;
    double energyPotential;

    bool isBound;
    bool isEjected;
};

#endif // PARTICLE_H
