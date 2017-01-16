#include <iostream>
#include <fstream>

#include "src/clustersystem.h"
#include "src/octree.h"
#include "src/barneshut.h"
#include "src/particle.h"

using namespace std;

int main()
{
    cout << "Hello Cluster!" << endl;

    ClusterSystem cs(1000, 40);

    cout << "Number of particles: " << cs.getNumberOfParticles() << endl << endl;

    cs.simulateBarnesHut(5, 0.001, true, "BH1000");
    //cs.simulateNtoN(5, 0.001, true, "NtoN1000");
    //cs.simulateNtoNRK4(5, 0.001, true, "NtoNRK1000");
    //cs.simulateNtoN(5, 0.001, true, "compareInt.smooth.NtoN100");
    //cs.simulateNtoNRK4(5, 0.001, true, "compareInt.smooth.NtoNRK100");

    return 0;
}

