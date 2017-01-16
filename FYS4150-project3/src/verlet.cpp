#include "verlet.h"

Verlet::Verlet(double x, double dx, double dt){
    /*
     * Constructor for Verlet integrator class
     *
     * @param double x, initial value
     * @param double dx, first derivative value
     * @param double dt, step length
     */
    this->iterationSteps = 0;

    this->stepLength = dt;
    this->xPrevPrev = x;
    // first step estimated with euler forward
    this->xPrev = this->xPrevPrev + dx * dt;

    this->iterationSteps++;
}

double Verlet::iterate(double ddx){
    /*
     * Advance one step forward
     *
     * @param double ddx, value for second derivative
     * @return double x, new value at time + dt
     */
    // first iteration uses euler
    if (this->iterationSteps == 1){
        this->iterationSteps++;
        return this->xPrev;
    }
    double xNew;
    xNew = 2 * this->xPrev - this->xPrevPrev + ddx * this->stepLength * this->stepLength;
    this->xPrevPrev = this->xPrev;
    this->xPrev = xNew;
    this->iterationSteps++;
    return xNew;
}

double Verlet::getVelocity(){
    /*
     * Get the average velocity over step.
     *
     * @return double v, velocity or first derivative
     */
    return (this->xPrev - this->xPrevPrev) / this->stepLength;
}

int Verlet::getNumberOfIterations(){
    /*
     * Return the current number of iterations.
     *
     * @return int number of iterations
     */
    return this->iterationSteps;
}
