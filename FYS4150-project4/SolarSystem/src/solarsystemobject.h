#ifndef SOLARSYSTEMOBJECT_H
#define SOLARSYSTEMOBJECT_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class SolarSystemObject{

public:
    SolarSystemObject(string, double, vec, vec);

    string getName();
    double getMass();
    vec getCurrentPosition();
    vec getCurrentVelocity();
    vec getForceG();
    double getForceGx();
    double getForceGy();
    double getForceGz();
    vec getDistanceToObject(SolarSystemObject);
    bool getHasLockedPosition();
    void setHasLockedPosition(bool);

    void setNewPosition(vec);
    void setNewPosition(double, double, double);
    void setNewVelocity(vec);
    void setNewVelocity(double, double, double);
    void increaseForceG(vec);
    void resetForceG();
private:
    string name;
    double mass;
    bool hasLockedPosition;

    vec currentPosition;
    vec currentVelocity;
    vec forceG;
};

#endif // SOLARSYSTEMOBJECT_H
