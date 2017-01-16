#ifndef EIGENJACOBIROTATION_H
#define EIGENJACOBIROTATION_H

#include<iostream>
#include<time.h>
#include<armadillo>

using namespace std;
using namespace arma;

class EigenJacobiRotation {

public:
    EigenJacobiRotation(mat);
    mat getRotatedMatrix();
    mat getEigenVectorMatrix();
    vec getEigenValues();
    int getMatrixSizen();
    int getIterationsNeeded();
    double getCalculationTime();
    double getTolerance();
    bool getEigenValuesFound();
    void findEigenValues();
    vec getFirstEigenVector();



private:
    void findMaxNonDiagValuePosition();
    void findMatrixDimensions();
    void rotateMatrix();
    void sortEigenvaluesMatrix();
    void findFirstEigenVector();

    mat A;
    mat AEigenVector;
    vec AEigenValues;
    vec FirstEigenVector;
    int n;
    double toleranceEpsilon;
    int iterationsNeeded;
    double calculationTime;
    bool eigenValuesFound;
    bool saveToFile;

    double maxNonDiagValue;
    int maxNonDiagValueColumn;
    int maxNonDiagValueRow;

};

#endif // EIGENJACOBIROTATION_H
