#include <iostream>

#include <src/solarsystem.h>

using namespace std;

int main(){
    cout << "Creating the solar system." << endl;
    SolarSystem solarsystem = SolarSystem(2, "solarsystemAU.txt");
    //SolarSystem solarsystem = SolarSystem(2, "sunearthAU.txt");
    solarsystem.createCelestialBody("sun", false);
    solarsystem.createCelestialBody("mercury");
    solarsystem.createCelestialBody("venus");
    solarsystem.createCelestialBody("earth");
    solarsystem.createCelestialBody("moon");
    solarsystem.createCelestialBody("mars");
    solarsystem.createCelestialBody("jupiter");
    solarsystem.createCelestialBody("saturn");
    solarsystem.createCelestialBody("uranus");
    solarsystem.createCelestialBody("neptune");
    solarsystem.createCelestialBody("pluto");

    //solarsystem.simulateVerlet(250.0, 0.01, "solarsystem.Verlet.dt0.01.dat");
    solarsystem.simulateRK4(250.0, 0.01, "solarsystem.RK4.dt0.01.dat");
    //solarsystem.simulateRK4(500.0, 0.005, "sun.earth.RK4.dt0.005.dat");
    //solarsystem.simulateVerlet(500, 0.01, "sun.earth.Verlet.dt0.01.dat");

    //SolarSystem solarsystem = SolarSystem(2, "sunearthjupiter1000AU.txt");
    //solarsystem.createCelestialBody("sun", false);
    //solarsystem.createCelestialBody("earth");
    //solarsystem.createCelestialBody("jupiter");
    //solarsystem.simulateVerlet(15.0, 0.001, "sun.earth.jupiter1000.Verlet.dt0.0001.unlocked.dat");
    //solarsystem.simulateRK4(15.0, 0.01, "sun.earth.jupiter1000.RK4.dt0.01.unlocked.dat");

    return 0;
}

