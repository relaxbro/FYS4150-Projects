#include "eigenjacobirotation.h"


EigenJacobiRotation::EigenJacobiRotation(mat matToBeRotated){
    /*
     * Constructor for class
     *
     * @param mat matToBeRotated
     */

    // set tolerance
    this->toleranceEpsilon = 1.0e-10;
    this->eigenValuesFound = false;

    this->A = matToBeRotated;
    // find matrix size
    findMatrixDimensions();
    // set up eigenvector matrix
    this->AEigenVector = eye<mat>(this->n, this->n);



}

vec EigenJacobiRotation::getEigenValues(){
    return this->AEigenValues;
}

mat EigenJacobiRotation::getRotatedMatrix(){
    return this->A;
}

mat EigenJacobiRotation::getEigenVectorMatrix(){
    return this->AEigenVector;
}

int EigenJacobiRotation::getMatrixSizen(){
    return this->n;
}

int EigenJacobiRotation::getIterationsNeeded(){
    return this->iterationsNeeded;
}

double EigenJacobiRotation::getCalculationTime(){
    return this->calculationTime;
}

double EigenJacobiRotation::getTolerance(){
    return this->toleranceEpsilon;
}

bool EigenJacobiRotation::getEigenValuesFound(){
    return this->eigenValuesFound;
}

vec EigenJacobiRotation::getFirstEigenVector(){
    return this->FirstEigenVector;
}

void EigenJacobiRotation::findFirstEigenVector(){
    /*
     * Find first eigenvector
     * Replaced by sortEigenvaluesMatrix()
     */

    int index = 0;
    for (int i = 0; i < this->n; i++){
        if (this->AEigenValues(0) == this->A(i, i)){
            index = i;
            break;
        }
    }
    this->FirstEigenVector = this->AEigenVector.col(index);
    return;
}

void EigenJacobiRotation::findEigenValues(){
    /*
     * Finding the eigenvalues
     */

    // check if allready found
    if (!(this->eigenValuesFound)){
        this->iterationsNeeded = 0;

        findMaxNonDiagValuePosition();

        // start timer
        clock_t timeStart, timeStop;
        timeStart = clock();

        // run until largest value is smaller than the tolerance
        while (this->maxNonDiagValue > this->toleranceEpsilon){

            findMaxNonDiagValuePosition();
            rotateMatrix();
            this->iterationsNeeded++;

        }
        // stop timer
        timeStop = clock();
        this->calculationTime = (timeStop - timeStart) / float(CLOCKS_PER_SEC);
        cout << "Calculation time: " << this->calculationTime << "s" << endl;
        cout << "Iterations needed: " << this->iterationsNeeded << endl;
        this->eigenValuesFound = true;

        // sort matrix
        sortEigenvaluesMatrix();
    }
    return;
}

void EigenJacobiRotation::rotateMatrix(){
    /*
     * Rotate the matrix
     */

    double sineValue = 0.0;
    double cosineValue = 1.0;

    // position of max value
    int k = this->maxNonDiagValueRow;
    int l = this->maxNonDiagValueColumn;

    if (this->A(k, l) != 0.0){
        // finding t and tau
        double tValue;
        double tauValue = (this->A(l, l) - this->A(k, k)) / (2 * this->A(k, l));

        if (tauValue < 0){
            tValue = 1.0 / (tauValue - sqrt(1.0 + tauValue * tauValue));
        }
        else{
            tValue = 1.0 / (tauValue + sqrt(1.0 + tauValue * tauValue));
        }
    // finding sine and cosine values
        cosineValue = 1.0 / sqrt(1.0 + tValue * tValue);
        sineValue = cosineValue * tValue;
    }

    double Akk = this->A(k, k);
    double All = this->A(l, l);
    double Akl = this->A(k, l);

    // change values at k, l
    this->A(k, k) = cosineValue * cosineValue * Akk - 2.0 * cosineValue * sineValue * Akl + sineValue * sineValue * All;
    this->A(l, l) = sineValue * sineValue * Akk + 2.0 * cosineValue * sineValue * Akl + cosineValue * cosineValue * All;
    this->A(k, l) = 0.0;
    this->A(l, k) = 0.0;

    double Aik, Ail;

    // change the rest of the matrix
    for (int i = 0; i < this->n; i++){
        if (i != k && i != l){
            Aik = this->A(i, k);
            Ail = this->A(i, l);
            this->A(i, k) = cosineValue * Aik - sineValue * Ail;
            this->A(i, l) = cosineValue * Ail + sineValue * Aik;
            this->A(k, i) = this->A(i, k);
            this->A(l, i) = this->A(i, l);
        }
        // updating the eigenvector matrix
        double AEigenVectorik = this->AEigenVector(i, k);
        double AEigenVectoril = this->AEigenVector(i, l);
        this->AEigenVector(i, k) = cosineValue * AEigenVectorik - sineValue * AEigenVectoril;
        this->AEigenVector(i, l) = cosineValue * AEigenVectoril + sineValue * AEigenVectorik;
    }
    return;
}

void EigenJacobiRotation::findMaxNonDiagValuePosition(){
    /*
     * Find max non-diagonal value and position
     */

    double maxValue = 0.0;
    int row = 0;
    int column = 0;

    // going thru the matrix
    for (int i = 0; i < n; i++){
        // only looking at upper part because of symmetry. Halves run time
        for (int j = i+1; j < n; j++){
            if (!(i == j)){
                if (abs(this->A(i, j)) > maxValue){
                    maxValue = abs(this->A(i, j));
                    row = i;
                    column = j;
                }
            }
        }
    }

    // setting value and position
    this->maxNonDiagValue = maxValue;
    this->maxNonDiagValueRow = row;
    this->maxNonDiagValueColumn = column;

    return;
}

void EigenJacobiRotation::findMatrixDimensions(){
    /*
     * Find dimension of matrix
     */
    int colums, rows;

    colums = this->A.n_cols;
    rows = this->A.n_rows;

    // must be a square matrix
    if (rows == colums){
        this->n = colums;
    }
    else{
        cout << "Size of matrix not n x n. Something went wrong." << endl;
        exit(0);
    }

    return;
}

void EigenJacobiRotation::sortEigenvaluesMatrix(){
    /*
     * Sort eigenvalues
     */

    // picking out eigenvalues and sort by value
    this->AEigenValues = sort(A.diag());

    int currentIndex = 0;
    int size = this->n;

    // going thru all eigenvalues and sorting matrix A and matrix AEigenVector
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            if (this->AEigenValues(i) == this->A(j, j)){
                currentIndex = j;
                break;
            }
        }
        // swap colum
        this->AEigenVector.swap_cols(i, currentIndex);
        this->A.swap_cols(i, currentIndex);
    }

    return;
}

