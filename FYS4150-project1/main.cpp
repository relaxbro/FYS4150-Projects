#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <armadillo>
#include <math.h>
#include <time.h>

using namespace std;
using namespace arma;


double sourceFunction(double x){
    /*
     * Returns the value of the source function f(x).
     *
     * @param double x, value to calculate f(x)
     * @return double function value f(x)
     */

    double f = 100 * exp(-10 * x);
    return f;
}

double ddSourceFunction(double x){
    /*
     * Returns the value of the double derivative f''(x)
     *
     * @param double x, value to calculate f''(x)
     * @return double function value f''(x)
     */

    double ddf = 1 - (1 - exp(-10)) * x - exp(-10 * x);
    return ddf;
}

void writeFile(int n, double h, vec v, vec relError, double relErrorMax, vec vLU = zeros<vec>(0)){
    /*
     * Writes data to file.
     * If no vLU is sent, create a vector of length 1. This to avoid sending vec if not wanted. Problems if set to NULL
     *
     * @param int n, steps (size of grid)
     * @param double h, step size
     * @param vec v, calculated values to be written
     * @param vec relError, relative errors to be ritten
     * @param vec vLU, solution from LU decomposition.
     * @return none
     */


    vec vLUnew;
    // if vLU is sent as parameter create a new vec with 0 on either end
    if (vLU.n_elem != 0){
        vLUnew.zeros(n+2);
        for (int i = 1; i < n+1; i++){
            vLUnew[i] = vLU[i-1];
        }
    }


    stringstream s;
    string filename;

    s << "data.n" << n << ".dat";

    filename = s.str();

    ofstream outfile;
    outfile.open(filename.c_str());

    // write steps n and step length h to file
    outfile << n << " " << h << " " << relErrorMax << endl;

    // write data to file
    if (vLU.n_elem == 0){
        for (int i = 0; i < n+2; i++){
            outfile << v[i] << " " << relError[i] << endl;
        }
    }
    else{
        for (int i = 0; i < n+2; i++){
            outfile << v[i] << " " << relError[i] <<  " " << vLUnew[i] << endl;
        }
    }

    // close the file
    outfile.close();
    cout << "File written: " << filename << "\n" << endl;
}

int main(int argc, char *argv[]){

    // missing command line arguemnts
    if (argc < 2){
        cout << "Command line arguments not correct!\nRun as: program integer -LU(if LU decomposition is wanted)" << endl;
        return 1;
    }

    // convert first argument(n) to int
    istringstream iss(argv[1]);
    int n;

    if (iss >> n){
       cout << "n = " << n << endl;
    }
    else{
        cout << "Argument not an integer." << endl;
        return 1;
    }

    clock_t tStart, tStop;

    // conditions and limits
    double lowerLimit = 0.0;
    double upperLimit = 1.0;
    double lowerBoundary = 0.0;
    double upperBoundary = 0.0;

    // step length
    double h = (upperLimit - lowerLimit) / (n + 1);
    cout << "Step length h = " << h << endl;

    // values in the tridiagonal metrix
    double a = -1;
    double b = 2;
    double c = -1;

    // vector for holding solution
    vec v(n+2);
    v.zeros();
    // vector for temporary values
    vec cPrime(n);
    cPrime.zeros();

    // setting initial values
    v[0] = lowerBoundary;
    v[n+1] = upperBoundary;
    // setting first "substitutor"
    double bTilde = b;
    v[1] = sourceFunction(h) / bTilde;

    // start timer for computational time
    tStart = clock();
    // substituting forwards
    for (int i = 2; i < n+1; i++){
        cPrime[i] = c / bTilde;
        bTilde = b - a * cPrime[i];
        v[i] = (sourceFunction(i * h) - a * v[i-1]) / bTilde;
    }
    // substituting backwards
    for (int i = n-1; i > 0; i--){
        v[i] = v[i] - cPrime[i+1] * v[i+1];
    }

    v = v * h * h;

    // stop timer
    tStop = clock();

    vec xValues = linspace<vec>(lowerLimit, upperLimit, n+2);
    // closed form solution
    vec uClosedForm;
    uClosedForm.zeros(n+2);

    vec relError;
    relError.zeros(n+2);

    // calculate the relative error
    for (int i = 1; i < n+1; i++){
        uClosedForm[i] = ddSourceFunction(xValues[i]);
        relError[i] = log10(abs((v[i] - uClosedForm[i]) / uClosedForm[i]));
    }

    double relErrorMax = relError[1];

    // find maximum relative error
    for (int i = 2; i < n+1; i++){
        if (relError[i] > relErrorMax){
            relErrorMax = relError[i];
        }
    }

    cout << "Max relative error: " << relErrorMax << endl;
    cout << "Calculation time: " << ((tStop - tStart) / float(CLOCKS_PER_SEC)) << "s" << endl;


    // check if user want LU decomposition
    if (argv[2]){
        if (string(argv[2]) != "-LU") {
            writeFile(n, h, v, relError, relErrorMax);
            return 0;
        }
    }
    // LU decomposition not wanted, finish program
    else{
        writeFile(n, h, v, relError, relErrorMax);
        return 0;
    }
    /*-------------------------------------------------------------
     *-------------------LU decomposition--------------------------
     *-------------------------------------------------------------
     */

    mat A, L, U;
    // try to allocate memory for matrices
    cout << "\nTrying to allocate memory for LU. If program abruptly shuts down no exception was thrown. Therefor no exception could be catched."  << endl;
    cout << "Remove -LU flag to not use LU decomposition." << endl;
    try{
        A.zeros(n, n);
        L.zeros(n, n);
        U.zeros(n, n);
    }
    // not enough memory, writing data so far
    // armadillo should throw std::bad_alloc (but does not allways)
    catch (std::bad_alloc){
        cout << "Not enough memory (RAM) for LU decomposition.\n" << endl;
        writeFile(n, h, v, relError, relErrorMax);
        return 1;
    }

    // generating the tridiagonal matrix A
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (j == i){
                A(i,j) = 2;
            }
            else if (j-i == 1){
                A(i,j) = -1;
            }
            else if (j-i == -1){
                A(i,j) = -1;
            }
        }
    }

    // values for source function
    vec source;
    source.zeros(n);
    for (int i = 0; i < n; i++){
        source[i] = sourceFunction((i+1)*h);
    }

    // start timer
    tStart = clock();

    // finding LU decomposition of A using armadillo lu
    lu(L, U, A);

    vec vLUtemp;
    vec vLU;

    // finding solution with armadillo solve
    vLUtemp = solve(L, source);
    vLU = solve(U, vLUtemp);
    vLU = vLU * h * h;

    // stop timer
    tStop = clock();

    cout << "\nCalculation time LU: " << ((tStop - tStart) / float(CLOCKS_PER_SEC)) << "s" << endl;

    writeFile(n, h, v, relError, relErrorMax, vLU);

    return 0;
}


