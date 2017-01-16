#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#include <armadillo>

using namespace arma;

class RungeKutta4{

public:
    RungeKutta4(vec, double);

    vec iterate(double);
private:
    int numberOfIterations;
    int numberOfCoupledDEs;

    double stepLength;
    vec prev;
    vec k1;
    vec k2;
    vec k3;
    vec k4;

};

#endif // RUNGEKUTTA4_H
