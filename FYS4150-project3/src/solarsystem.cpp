#include "solarsystem.h"


SolarSystem::SolarSystem(int dimension, string bodiesFilename){
    /*
     * Construcot for SolarSystem class
     *
     * @param int dimension
     * @param string bodiesFilename, filename where initial values of objects/bodies are stored
     */
    this->dimension = dimension;
    this->massSystem = 0;
    this->currentAngulerMomentum = zeros<vec>(3);
    // set the gravitation constant
    this->gravitationalConstant = 4 * M_PI * M_PI;
    // file in initialdata folder
    this->bodiesFilename = "initialdata/" + bodiesFilename;
}

int SolarSystem::getDimension(){
    return this->dimension;
}

int SolarSystem::getNumberOfCelestialBodies(){
    return this->celestialBodies.size();
}

void SolarSystem::addCelestialBody(SolarSystemObject body){
    /*
     * Add new celestial body to list/vector.
     *
     * @param SolarSystemObject body
     */
    celestialBodies.push_back(body);
    return;
}

void SolarSystem::createCelestialBody(string name, bool locked){
    /*
     * Create new celestial body.
     *
     * @param string name
     * @param bool locked, if position is locked
     */
    // open data file with initial values
    ifstream infile;
    try {
        infile.open(this->bodiesFilename.c_str());
    }
    catch(...){
        cout << "Failed to open file: " << this->bodiesFilename << endl;
        return;
    }

    double posx, posy, posz, velx, vely, velz, mass;
    bool bodyFound = false;
    string line;

    // read data from file and convert
    while (getline(infile, line)){
        vector<string> words;
        // split line at whitespace
        boost::split(words, line, boost::is_any_of(" "));
        if (words[0] == name){
            posx = atof(words[1].c_str());
            posy = atof(words[2].c_str());
            posz = atof(words[3].c_str());
            velx = atof(words[4].c_str());
            vely = atof(words[5].c_str());
            velz = atof(words[6].c_str());
            mass = atof(words[7].c_str());
            bodyFound = true;
            break;
        }
    }

    // celestial body not found in file
    if(!bodyFound){
        cout << name << " not found in datafile: " << this->bodiesFilename << " . Make sure it is there." << endl;
        return;
    }

    vec pos = zeros<vec>(3);
    pos(0) = posx;
    pos(1) = posy;
    pos(2) = posz;
    vec vel = zeros<vec>(3);
    vel(0) = velx;
    vel(1) = vely;
    vel(2) = velz;
    cout << "--- Adding" << endl;
    cout << name << ", mass: " << mass << endl;
    cout << "position: " << posx << " " << posy << " " << posz << endl;
    // create object/body
    SolarSystemObject body(name, mass, pos, vel);

    if (locked){
        body.setHasLockedPosition(locked);
    }

    addCelestialBody(body);
    // adding mass to total mass
    this->massSystem += mass;
    return;
}

void SolarSystem::findForceOnBody(int indexOfBody){
    /*
     * Find all forces on body i.
     *
     * @param int indexOfBody, index of body in list
     */
    // pointer to object
    SolarSystemObject &currentObject = this->celestialBodies[indexOfBody];
    // find all forces
    for (int j = 0; j < this->getNumberOfCelestialBodies(); j++){
        if (j == indexOfBody){
            continue;
        }
        vec r = currentObject.getDistanceToObject(this->celestialBodies[j]);

        double GMM = this->gravitationalConstant * currentObject.getMass() * this->celestialBodies[j].getMass();
        vec force = - r * GMM / pow(norm(r), 3);

        currentObject.increaseForceG(force);
    }
    return;
}

void SolarSystem::findForceOnAllBodies(){
    /*
     * Find all forces working on all bodies
     */
    this->resetForceOnAllBodies();
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        for (int j = i+1; j < this->getNumberOfCelestialBodies(); j++){
            vec r = currentObject.getDistanceToObject(this->celestialBodies[j]);

            double GMM = this->gravitationalConstant * currentObject.getMass() * this->celestialBodies[j].getMass();
            vec force = - r * GMM / pow(norm(r), 3);

            currentObject.increaseForceG(force);
            // add corresponding force to other object
            this->celestialBodies[j].increaseForceG(-force);
        }
    }
    return;
}

void SolarSystem::resetForceOnAllBodies(){
    /*
     * Reset all forces back to zero.
     */
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        currentObject.resetForceG();
    }
    return;
}

void SolarSystem::findEnergyKinetic(){
    /*
     * Find the total kinetic energy of the system.
     */
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        double vel2 = pow(norm(currentObject.getCurrentVelocity()), 2);
        this->currentEnergyKinetic += 0.5 * currentObject.getMass() * vel2;
    }
    return;
}

void SolarSystem::findEnergyPotential(){
    /*
     * Find the total potential energy of the system.
     */
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObjet = this->celestialBodies[i];
        for (int j = 0; j < this->getNumberOfCelestialBodies(); j++){
            if (j == i){
                continue;
            }
            SolarSystemObject &tempObject = this->celestialBodies[j];
            double r = norm(currentObjet.getDistanceToObject(tempObject));
            this->currentEnergyPotential += - this->gravitationalConstant * currentObjet.getMass() * tempObject.getMass() / r;
        }
    }
    /*
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        double r = norm(currentObject.getCurrentPosition());
        if (r > 1e-5){
            this->currentEnergyPotential -= this->gravitationalConstant * this->massSystem * currentObject.getMass() / r;
        }
    }
    */
    return;
}

void SolarSystem::resetEnergy(){
    /*
     * Reset kinetic and potential energy to zero.
     */
    this->currentEnergyKinetic = 0;
    this->currentEnergyPotential = 0;
    return;
}

void SolarSystem::findAngularMomentum(){
    /*
     * Find the angular momentum of the system.
     */
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        this->currentAngulerMomentum += cross(currentObject.getCurrentPosition(), currentObject.getMass() * currentObject.getCurrentVelocity());
    }
    return;
}

void SolarSystem::resetAngularMomentum(){
    /*
     * Reset angular momentum to zero.
     */
    this->currentAngulerMomentum(0) = 0;
    this->currentAngulerMomentum(1) = 0;
    this->currentAngulerMomentum(2) = 0;
    return;
}

void SolarSystem::simulateVerlet(double t, double dt, string outputFilename){
    /*
     * Simulate for a time t using Verlet.
     *
     * @param double t, simulation time
     * @param double dt, step length
     * @param string outputFilename, base for naming all outputfiles
     */
    this->simulationTime = t;
    this->simulationTimeCurrent = 0;
    this->simulationTimeStep = dt;
    this->simulationSteps = this->simulationTime / this->simulationTimeStep;
    this->numberOfIntegrationSteps = 0;

    // start timer
    clock_t timeStart, timeStop;
    timeStart = clock();

    // list of integrators
    vector<Verlet> xIntegration;
    vector<Verlet> yIntegration;
    vector<Verlet> zIntegration;
    // initiate the integrators
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
       SolarSystemObject &currentObject = this->celestialBodies[i];
        vec pos = currentObject.getCurrentPosition();
        vec vel = currentObject.getCurrentVelocity();
        double x = pos(0);
        double y = pos(1);
        double z = pos(2);
        double dx = vel(0);
        double dy = vel(1);
        double dz = vel(2);
        Verlet tempx = Verlet(x, dx, this->simulationTimeStep);
        Verlet tempy = Verlet(y, dy, this->simulationTimeStep);
        Verlet tempz = Verlet(z, dz, this->simulationTimeStep);
        xIntegration.push_back(tempx);
        yIntegration.push_back(tempy);
        zIntegration.push_back(tempz);
    }

    // finding total kinetic and potential energy
    this->resetEnergy();
    this->findEnergyKinetic();
    this->findEnergyPotential();

    // find total angular momentum
    this->resetAngularMomentum();
    this->findAngularMomentum();

    // open outputfiles and write initial data
    this->openOutputFiles(outputFilename);


    //------- integration ---------
    cout << "\nStarting integration Verlet" << endl;
    while (this->simulationTimeCurrent < this->simulationTime){
        this->simulationTimeCurrent += this->simulationTimeStep;

        // find forces
        this->resetForceOnAllBodies();
        this->findForceOnAllBodies();

        // advance one step for each object
        for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
            SolarSystemObject &currentObject = this->celestialBodies[i];
            vec force = currentObject.getForceG();
            double mass = currentObject.getMass();
            double ddx = force(0) / mass;
            double ddy = force(1) / mass;
            double ddz = force(2) / mass;
            double x = xIntegration[i].iterate(ddx);
            double y = yIntegration[i].iterate(ddy);
            double z = zIntegration[i].iterate(ddz);
            double dx = xIntegration[i].getVelocity();
            double dy = yIntegration[i].getVelocity();
            double dz = zIntegration[i].getVelocity();

            // update velocit and position
            currentObject.setNewPosition(x, y, z);
            currentObject.setNewVelocity(dx, dy, dz);

        }
        // find energy
        this->resetEnergy();
        this->findEnergyKinetic();
        this->findEnergyPotential();

        // find total angular momentum
        this->resetAngularMomentum();
        this->findAngularMomentum();

        // write to file
        this->writeOutputFiles();
    }
    // close file
    this->closeOutputFiles();

    // stop timer
    timeStop = clock();
    cout << "Simulation time: " << ((timeStop - timeStart) / float(CLOCKS_PER_SEC)) << "s." << endl;
    return;
}

void SolarSystem::simulateRK4(double t, double dt, string outputFilename){
    /*
     * Simulate for a time t using Runge-Kutta 4.
     *
     * @param double t, simulation time
     * @param double dt, step length
     * @param string outputFilename, base for naming all outputfiles
     */
    this->simulationTime = t;
    this->simulationTimeCurrent = 0;
    this->simulationTimeStep = dt;
    this->simulationSteps = this->simulationTime / this->simulationTimeStep;
    this->numberOfIntegrationSteps = 0;

    // start timer
    clock_t timeStart, timeStop;
    timeStart = clock();

    // finding total kinetic and potential energy
    this->resetEnergy();
    this->findEnergyKinetic();
    this->findEnergyPotential();

    // find total angular momentum
    this->resetAngularMomentum();
    this->findAngularMomentum();

    // open outputfiles and write initial data
    this->openOutputFiles(outputFilename);

    //------- integration ---------
    cout << "\nStarting integration RK4" << endl;
    // here goes while loop for integration
    while (this->simulationTimeCurrent < this->simulationTime){
        this->simulationTimeCurrent += this->simulationTimeStep;

        // find all forces
        this->findForceOnAllBodies();

        // advance one step
        this->RK4iterate(this->simulationTimeStep);

        // find energy
        this->resetEnergy();
        this->findEnergyKinetic();
        this->findEnergyPotential();

        // find total angular momentum
        this->resetAngularMomentum();
        this->findAngularMomentum();

        // write to file
        this->writeOutputFiles();
    }
    // close file
    this->closeOutputFiles();

    // stop timer
    timeStop = clock();
    cout << "Simulation time: " << ((timeStop - timeStart) / float(CLOCKS_PER_SEC)) << "s." << endl;
    return;
}

void SolarSystem::RK4iterate(double dt){
    /*
     * Advance one step with Runge-Kutta 4.
     *
     * @param double dt, step size/length
     */
    int n = this->getNumberOfCelestialBodies();

    // matrices for storing intermediate step values
    mat k1 = zeros<mat>(3, n);
    mat k2 = zeros<mat>(3, n);
    mat k3 = zeros<mat>(3, n);
    mat k4 = zeros<mat>(3, n);

    mat l1 = zeros<mat>(3, n);
    mat l2 = zeros<mat>(3, n);
    mat l3 = zeros<mat>(3, n);
    mat l4 = zeros<mat>(3, n);

    // storing velocity and position at current time
    mat currentVelocity = zeros<mat>(3, n);
    mat currentPosition = zeros<mat>(3, n);

    // calculate l1 and k1
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        currentPosition.col(i) = currentObject.getCurrentPosition();
        currentVelocity.col(i) = currentObject.getCurrentVelocity();

        l1.col(i) = dt * currentObject.getForceG() / currentObject.getMass();
        k1.col(i) = dt * currentObject.getCurrentVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l1.col(i) / 2.0;
        vec newPos = currentPosition.col(i) + k1.col(i) / 2.0;
        currentObject.setNewVelocity(newVel);
        currentObject.setNewPosition(newPos);
    }

    this->findForceOnAllBodies();

    // calculate l2 and k2
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        l2.col(i) = dt * currentObject.getForceG() / currentObject.getMass();
        k2.col(i) = dt * currentObject.getCurrentVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l2.col(i) / 2.0;
        vec newPos = currentPosition.col(i) + k2.col(i) / 2.0;
        currentObject.setNewVelocity(newVel);
        currentObject.setNewPosition(newPos);
    }

    this->findForceOnAllBodies();

    // calculate l3 and l4
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        l3.col(i) = dt * currentObject.getForceG() / currentObject.getMass();
        k3.col(i) = dt * currentObject.getCurrentVelocity();

        // set new temp velocity and position
        vec newVel = currentVelocity.col(i) + l3.col(i);
        vec newPos = currentPosition.col(i) + k3.col(i);
        currentObject.setNewVelocity(newVel);
        currentObject.setNewPosition(newPos);

    }

    this->findForceOnAllBodies();

    // calculate l4 and k4
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        l4.col(i) = dt * currentObject.getForceG() / currentObject.getMass();
        k4.col(i) = dt * currentObject.getCurrentVelocity();

        // set final new velocity and position
        vec newVel = currentVelocity.col(i) + (l1.col(i) + 2*l2.col(i) + 2*l3.col(i) + l4.col(i)) / 6.0;
        vec newPos = currentPosition.col(i) + (k1.col(i) + 2*k2.col(i) + 2*k3.col(i) + k4.col(i)) / 6.0;
        currentObject.setNewVelocity(newVel);
        currentObject.setNewPosition(newPos);
    }
    return;
}

void SolarSystem::openOutputFiles(string outputFilename){
    /*
     * Open all outfiles and write initial data.
     *
     * @param string outputFilename, base for filenames for all output files
     */
    string outputFilenamePosition = "position." + outputFilename;
    string outputFilenameVelocity = "velocity." + outputFilename;
    string outputFilenameForce = "force." + outputFilename;
    string outputFilenameEnergy = "energy." + outputFilename;
    string outputFilenameAngularMomentum = "angular.momentum." + outputFilename;
    // open files
    outfilePosition = new ofstream(outputFilenamePosition.c_str());
    outfileVelocity = new ofstream(outputFilenameVelocity.c_str());
    outfileForce = new ofstream(outputFilenameForce.c_str());
    outfileEnergy = new ofstream(outputFilenameEnergy.c_str());
    outfileAngularMomentum = new ofstream(outputFilenameAngularMomentum.c_str());

    // write initial data for position
    *outfilePosition << "Number of celestial bodies: " << this->getNumberOfCelestialBodies() << endl;
    *outfilePosition << "Time: " << this->simulationTime << endl;
    *outfilePosition << "dt: " << this->simulationTimeStep << endl;
    *outfilePosition << "Objects(mass): ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        *outfilePosition << this->celestialBodies[i].getName() << " " << this->celestialBodies[i].getMass() << " ";
    }
    *outfilePosition << endl;

    *outfilePosition << "time -- x1 -- y1 -- z1 -- x2 -- y2 -- z2 -- ..." << endl;
    *outfilePosition << 0 << " ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        vec pos = this->celestialBodies[i].getCurrentPosition();
        *outfilePosition << pos(0) << " " << pos(1) << " " << pos(2) << " ";
    }
    *outfilePosition << endl;

    // write initial data for velocity
    *outfileVelocity << "Number of celestial bodies: " << this->getNumberOfCelestialBodies() << endl;
    *outfileVelocity << "Time: " << this->simulationTime << endl;
    *outfileVelocity << "dt: " << this->simulationTimeStep << endl;
    *outfileVelocity << "Objects: ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        *outfileVelocity << this->celestialBodies[i].getName() << " -- ";
    }
    *outfileVelocity << endl;

    *outfileVelocity << "time -- dx1 -- dy1 -- dz1 -- dx2 -- dy2 -- dz2 -- ..." << endl;
    *outfileVelocity << 0 << " ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        vec vel = this->celestialBodies[i].getCurrentVelocity();
        *outfileVelocity << vel(0) << " " << vel(1) << " " << vel(2) << " ";
    }
    *outfileVelocity << endl;

    // write initial data for force
    *outfileForce << "Number of celestial bodies: " << this->getNumberOfCelestialBodies() << endl;
    *outfileForce << "Time: " << this->simulationTime << endl;
    *outfileForce << "dt: " << this->simulationTimeStep << endl;
    *outfileForce << "Objects: ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        *outfileForce << this->celestialBodies[i].getName() << " -- ";
    }
    *outfileForce << endl;

    *outfileForce << "time -- Fx1 -- Fy1 -- Fz1 -- Fx2 -- Fy2 -- Fz2 -- ..." << endl;
    *outfileForce << 0 << " ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        vec force = this->celestialBodies[i].getForceG();
        *outfileForce << force(0) << " " << force(1) << " " << force(2) << " ";
    }
    *outfileForce << endl;

    // write initial data for energy
    *outfileEnergy << "Number of celestial bodies: " << this->getNumberOfCelestialBodies() << endl;
    *outfileEnergy << "Time: " << this->simulationTime << endl;
    *outfileEnergy << "dt: " << this->simulationTimeStep << endl;
    *outfileEnergy << "Objects: ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        *outfileEnergy << this->celestialBodies[i].getName() << " -- ";
    }
    *outfileEnergy << endl;

    *outfileEnergy << "time -- kinetic -- potential -- total" << endl;
    *outfileEnergy << 0 << " ";
    *outfileEnergy << this->currentEnergyKinetic << " " << this->currentEnergyPotential << " " << (this->currentEnergyKinetic + this->currentEnergyPotential);
    *outfileEnergy << endl;

    // write initial data for angular momentum
    *outfileAngularMomentum << "Number of celestial bodies: " << this->getNumberOfCelestialBodies() << endl;
    *outfileAngularMomentum << "Time: " << this->simulationTime << endl;
    *outfileAngularMomentum << "dt: " << this->simulationTimeStep << endl;
    *outfileAngularMomentum << "Objects: ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        *outfileAngularMomentum << this->celestialBodies[i].getName() << " -- ";
    }
    *outfileAngularMomentum << endl;

    *outfileAngularMomentum << "time -- Lx -- Ly -- Lz" << endl;
    *outfileAngularMomentum << 0 << " ";
    *outfileAngularMomentum << this->currentAngulerMomentum(0) << " " << this->currentAngulerMomentum(1) << " " << this->currentAngulerMomentum(2);
    *outfileAngularMomentum << endl;

    return;
}

void SolarSystem::writeOutputFiles(){
    /*
     * Write data to file.
     */
    *outfilePosition << this->simulationTimeCurrent << " ";
    *outfileVelocity << this->simulationTimeCurrent << " ";
    *outfileForce << this->simulationTimeCurrent << " ";
    *outfileEnergy << this->simulationTimeCurrent << " ";
    *outfileAngularMomentum << this->simulationTimeCurrent << " ";
    for (int i = 0; i < this->getNumberOfCelestialBodies(); i++){
        SolarSystemObject &currentObject = this->celestialBodies[i];
        *outfilePosition << currentObject.getCurrentPosition()(0) << " " << currentObject.getCurrentPosition()(1) << " " << currentObject.getCurrentPosition()(2) << " ";
        *outfileVelocity << currentObject.getCurrentVelocity()(0) << " " << currentObject.getCurrentVelocity()(1) << " " << currentObject.getCurrentVelocity()(2) << " ";
        *outfileForce << currentObject.getForceG()(0) << " " << currentObject.getForceG()(1) << " " << currentObject.getForceG()(2) << " ";
    }
    *outfileEnergy << this->currentEnergyKinetic << " " << this->currentEnergyPotential << " " << (this->currentEnergyKinetic + this->currentEnergyPotential);
    *outfileAngularMomentum << this->currentAngulerMomentum(0) << " " << this->currentAngulerMomentum(1) << " " << this->currentAngulerMomentum(2);
    // add newline
    *outfilePosition << endl;
    *outfileVelocity << endl;
    *outfileForce << endl;
    *outfileEnergy << endl;
    *outfileAngularMomentum << endl;

    return;
}

void SolarSystem::closeOutputFiles(){
    /*
     * Close all outfiles.
     */
    outfilePosition->close();
    outfileVelocity->close();
    outfileForce->close();
    outfileEnergy->close();
    outfileAngularMomentum->close();

    return;
}
