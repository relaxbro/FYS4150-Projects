#ifndef VERLET_H
#define VERLET_H

#include<iostream>

class Verlet{
public:
    Verlet(double, double, double);

    double iterate(double);
    double getVelocity();
    int getNumberOfIterations();
private:
    int iterationSteps;
    double xPrev;
    double xPrevPrev;
    double stepLength;
};

#endif // VERLET_H
