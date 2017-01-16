#include "rungekutta4.h"

RungeKutta4::RungeKutta4(vec input, double stepLength){
    this-> prev = input;
    this->stepLength = stepLength;
}

vec RungeKutta4::iterate(double t){
    //k1 = stepLenght * f(prev, t);
    //k2 = stepLenght * f(prev + k1/2.0, t + this->stepLength/2);
    //k3 = stepLenght * f(prev + k2/2.0, t + this->stepLength/2);
    //k4 = stepLenght * f(prev + k3, t + this->steplength);
    //prev += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    return prev;
}
