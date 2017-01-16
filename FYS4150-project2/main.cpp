#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <math.h>
#include <time.h>
#include <vector>

#include "src/eigenjacobirotation.h"

using namespace std;
using namespace arma;

void schroedinger(int n, int numberOfElectrons, double rhoFinal, double omega);
void writeToFile(EigenJacobiRotation sup, int n, int numberOfElectrons, double rhoFinal, double omega);
void iterationsVSDimensionality(int nMin, int nMax, int dn, double rhoFinal);
void armadilloCompare(int nMin, int nMax, int dn, double rhoFinal);
void multipleOmega(int n, int numberOfElectrons, double rhoFinal);

int main(int argc, char *argv[]){

    // running for multiple values at once. used to create lots of data
    if (false){
        //iterationsVSDimensionality(10, 400, 10, 5.0);
        //armadilloCompare(10, 400, 10, 5.0);
        multipleOmega(200, 2, 5.0);
        return 0;
    }


    // missing command line arguments
    if (argc < 4){
        cout << "Missing arguments!\nRun as: program n(int) electrons(int) rhoFinal(double) omega(double)\n" << endl;
        return 0;
    }

    int n;
    int numberOfElectrons;
    double rhoFinal;
    double omega = 0;

    // convert first argument
    istringstream iss1(argv[1]);
    if (iss1 >> n){
        cout << "n = " << n;
    }
    else{
        cout << "First argument not integer" << endl;
        return 0;
    }

    // convert second argument
    istringstream iss2(argv[2]);
    if (iss2 >> numberOfElectrons){
        cout << ", numberOfElectrons = " << numberOfElectrons;
    }
    else{
        cout << "\nSecond argument not integer" << endl;
        return 0;
    }

    // convert third argument
    istringstream iss3(argv[3]);
    if (iss3 >> rhoFinal){
        cout << ", rhoFinal = " << rhoFinal;
    }
    else{
        cout << "\nThird argument not integer or double" << endl;
        return 0;
    }

    if (argc > 4){
        // convert fourth argument
        istringstream iss4(argv[4]);
        if (iss4 >> omega){
            cout << ", omega = " << omega << endl;
        }
        else{
            cout << "\nFourth argument not integer or double" << endl;
            return 0;
        }
    }

    if (numberOfElectrons != 1 && numberOfElectrons != 2){
        cout << "\nOnly 1 or 2 electrons possible. numberOfElectrons = " << numberOfElectrons << endl;
        return 0;
    }

    if (numberOfElectrons == 2 && omega == 0){
        cout << "Oscillator strength cannot be 0 with two electrons. omega = " << omega << endl;
    }

    // solving the problem
    schroedinger(n, numberOfElectrons, rhoFinal, omega);

    cout << endl;
    return 0;
}

void schroedinger(int n, int numberOfElectrons, double rhoFinal, double omega){
    /*
     * Setting up problem and solving with jacobi method
     *
     * @param int n, size of matrix
     * @param int numberOfElectrons
     * @param double rhoFinal, max value for rho
     * @param double omega, frquenzy
     * @return none
     */

    // setting inital rho
    double rhoStart = 0.0;
    double h = (rhoFinal - rhoStart) / n;

    // initializing matrix and vector
    mat A = zeros<mat>(n, n);
    vec V = zeros<vec>(n);
    double rho;

    // setting up potential for one electron
    if (numberOfElectrons == 1){
        for (int i = 0; i < n; i++){
            rho = rhoStart + (i + 1) * h;
            V(i) = rho * rho;
        }
    }
    // setting up potential for two electrons
    else {
        for (int i = 0; i < n; i++){
            rho = rhoStart + (i +1) * h;
            V(i) = (omega * omega) * (rho * rho) + 1.0 / rho;
        }
    }

    // vector with values over and under diagonal
    vec e = zeros<vec>(n-1);
    e += - 1.0 / (h * h);

    // inserting values in matrix
    A.diag(1) = e;
    A.diag() = 2.0 / (h * h) + V;
    A.diag(-1) = e;

    // creating object
    EigenJacobiRotation sup = EigenJacobiRotation(A);
    // find eigenvalues
    sup.findEigenValues();
    // write to file
    writeToFile(sup, n, numberOfElectrons, rhoFinal, omega);

    return;
}


void writeToFile(EigenJacobiRotation sup, int n, int numberOfElectrons, double rhoFinal, double omega){
    /*
     * Writing data to file
     *
     * @param EigenJacobiRotation sup
     * @param int n, size of matrix
     * @param int numberOfElectrons
     * @param double rhoFinal, max value for rho
     * @param double omega, frquenzy
     * @return none
     */

    // creating file name
    stringstream s;
    s << "data.e" << numberOfElectrons << ".n" << n << ".rho" << rhoFinal << ".w" << omega << ".dat";
    string filename = s.str();

    ofstream outfile;
    // open file
    outfile.open(filename.c_str());

    // write to file
    outfile << "Parameters (omega ignored when numberOfElectrons = 1):" << endl;
    outfile << "n: " << n << endl;
    outfile << "numberOfElectrons: " << numberOfElectrons << endl;
    outfile << "rhoFinal: " << rhoFinal << endl;
    outfile << "omega: " << omega << endl;
    outfile << "Tolerance: " << sup.getTolerance() << endl;
    outfile << "Iterations needed: " << sup.getIterationsNeeded() << endl;
    outfile << "Calculation time; " << sup.getCalculationTime() << endl << endl;

    outfile << "Eigenvalues--- ---First Eigenvector--- ---Second--- ---Third" << endl;

    for (int i = 0; i < n; i++){
        outfile << sup.getEigenValues()(i) << " " << sup.getEigenVectorMatrix().col(0)(i) << " ";
        outfile << sup.getEigenVectorMatrix().col(1)(i) << " " << sup.getEigenVectorMatrix().col(2)(i) << endl;
    }

    // close file
    outfile.close();
    cout << "File written: " << filename << endl;
    return;
}


void iterationsVSDimensionality(int nMin, int nMax, int dn, double rhoFinal){
    /*
     * Function to run for multiple matrix sizes n.
     *
     * @param int nMin, minimum n
     * @param int nMax, maximum n
     * @param int dn, stepsize for n
     * @param double rhoFinal, max value for rho
     * @return none
     */

    // number of steps
    int nSteps = (nMax - nMin) / dn;

    vec n = zeros<vec>(nSteps+1);

    for (int i = 0; i < nSteps+1; i++){
        n(i) = (i + 1) * dn;
    }

    double rhoStart = 0.0;

    // setting up vector containing multiple EigenJacobiRotation objects
    vector<EigenJacobiRotation> sup;
    // reserve for spped up
    sup.reserve(nSteps+1);
    // solve for all values of n
    for (int i = 0; i < nSteps+1; i++){

        int nTemp = n(i);
        cout << "n = " << nTemp << endl;
        double h = (rhoFinal - rhoStart) /nTemp;
        // creating matrix and vector
        mat A = zeros<mat>(nTemp, nTemp);
        vec V = zeros<vec>(nTemp);
        double rho;

        // find potential
        for (int j = 0; j < nTemp; j++){
            rho = rhoStart + (j + 1) * h;
            V(j) = rho * rho;
        }

        // vector with values over and under diagonal
        vec e = zeros<vec>(nTemp-1);
        e += - 1.0 / (h * h);

        // inserting values into matrix
        A.diag(1) = e;
        A.diag() = 2.0 / (h * h) + V;
        A.diag(-1) = e;

        // creating object in object list
        sup.push_back(EigenJacobiRotation(A));
        // solving
        sup[i].findEigenValues();

    }

    // creating file name
    stringstream s;
    s << "data.n" << nMin << "." << nMax << ".rhoFinal" << rhoFinal << ".dat";
    string filename = s.str();

    // open file
    ofstream outfile;
    outfile.open(filename.c_str());

    // write data to file
    outfile << "Comparing number of iterations needed." << endl;
    outfile << "nMin: " << nMin << endl;
    outfile << "nMax: " << nMax << endl;
    outfile << "dn: " << dn << endl;
    outfile << "rhoFianl: " << rhoFinal << endl;
    outfile << "Tolerance: " << sup[0].getTolerance() << endl << endl;

    outfile << "---n--- ---iterations--- ---calculationTime--- ---EV1--- ---EV2--- ---EV3--- ---EV4---" << endl;

    for (int i = 0; i < nSteps+1; i++){
        outfile << sup[i].getMatrixSizen() << " " << sup[i].getIterationsNeeded() << " ";
        outfile << sup[i].getCalculationTime() << " ";
        outfile << sup[i].getEigenValues()(0) << " " << sup[i].getEigenValues()(1) << " ";
        outfile << sup[i].getEigenValues()(2) << " " << sup[i].getEigenValues()(3) << endl;
    }

    // close file
    outfile.close();
    cout << "File written: " << filename << endl;

    return;
}

void armadilloCompare(int nMin, int nMax, int dn, double rhoFinal){
    /*
     * Function to solve problem with armadillo eig_sym
     *
     * @param int nMin, minimum n
     * @param int nMax, maximum n
     * @param int dn, stepsize for n
     * @param double rhoFinal, max value for rho
     * @return none
     */

    cout << "\n------------- Using armadillo" << endl;

    vec eigval;
    mat eigvec;
    clock_t timeStart, timeStop;
    // number of steps
    int nSteps = (nMax - nMin) / dn;

    vec n = zeros<vec>(nSteps+1);

    for (int i = 0; i < nSteps+1; i++){
        n(i) = (i + 1) * dn;
    }

    double rhoStart = 0.0;

    // create file name
    stringstream s;
    s << "data.n" << nMin << "." << nMax << ".rhoFinal" << rhoFinal << ".arma.dat";
    string filename = s.str();

    // open file
    ofstream outfile;
    outfile.open(filename.c_str());

    // write data to file
    outfile << "Using arma::eig_sym to find values" << endl;
    outfile << "nMin: " << nMin << endl;
    outfile << "nMax: " << nMax << endl;
    outfile << "dn: " << dn << endl;
    outfile << "rhoFinal: " << rhoFinal << endl << endl;

    outfile << "---n--- ---calculationTime--- ---EV1--- ---EV2--- ---EV3--- ---EV4---" << endl;

    // run for all n
    for (int i = 0; i < nSteps+1; i++){

        int nTemp = n(i);
        double h = (rhoFinal - rhoStart) /nTemp;

        // create matrix and vector
        mat A = zeros<mat>(nTemp, nTemp);
        vec V = zeros<vec>(nTemp);
        double rho;

        // find values for potential
        for (int j = 0; j < nTemp; j++){
            rho = rhoStart + (j + 1) * h;
            V(j) = rho * rho;
        }

        vec e = zeros<vec>(nTemp-1);
        e += - 1.0 / (h * h);

        // insert values into matrix
        A.diag(1) = e;
        A.diag() = 2.0 / (h * h) + V;
        A.diag(-1) = e;

        // start timer
        timeStart = clock();
        // solve with armadillo
        eig_sym(eigval, eigvec, A);
        // end timer
        timeStop = clock();
        double calculationTime = (timeStop - timeStart) / float(CLOCKS_PER_SEC);

        cout << "n: " << nTemp << " Calculation time: " << calculationTime << " s" << endl;

        // write wanted values to file
        outfile << nTemp << " " << calculationTime << " ";
        outfile << eigval(0) << " " << eigval(1) << " " << eigval(2) << " " << eigval(3) << endl;

    }

    // close file
    outfile.close();

    return;
}

void multipleOmega(int n, int numberOfElectrons, double rhoFinal){
    /*
     * Function to run for multiple vales of omega.
     *
     * @param int n
     * @param int numberOfElectrons
     * @param double rhoFinal
     * @return none
     */

    // number of steps for omega
    int omegaSteps = 10;
    // setting up vector for omerga
    vec omega = zeros<vec>(omegaSteps);
    omega(0) = 0.01;
    omega(1) = 0.05;
    omega(2) = 0.1;
    omega(3) = 0.25;
    omega(4) = 0.5;
    omega(5) = 1.0;
    omega(6) = 1.5;
    omega(7) = 2.5;
    omega(8) = 5.0;
    omega(9) = 10.0;

    // solve for all values of omega
    for (int i = 0; i < omegaSteps; i++){
        schroedinger(n, numberOfElectrons, rhoFinal, omega(i));
    }

    return;
}

