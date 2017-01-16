import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True})

class EnergyData:

    def __init__(self, N, letter):
        self.time = []
        self.Ep = []
        self.Ek = []
        self.Etot = []
        self.bound = []
        self.unbound = []
        self.ejected = []
        self.N = N
        self.filename = 'data/energyBH%d%s.dat' %(self.N, letter)
        try:  
            fileEnergy = open(self.filename, 'r')
        except: 
            print '%s not found.' %(self.filename)
            sys.exit()

        for line in fileEnergy:
            words = line.split()
            self.time.append(float(words[0]))
            self.Ek.append(float(words[1]))
            self.Ep.append(float(words[2]))
            self.bound.append(float(words[3]))
            self.ejected.append(float(words[4]))
        fileEnergy.close()

        self.time = np.asarray(self.time)
        self.Ek = np.asarray(self.Ek)
        self.Ep = np.asarray(self.Ep)
        self.Etot = self.Ek + self.Ep
        self.bound = np.asarray(self.bound)
        self.ejected = np.asarray(self.ejected)
        self.unbound = self.N - self.bound

energy100a = EnergyData(100, 'a')
energy100b = EnergyData(100, 'b')
energy100c = EnergyData(100, 'c')
energy100d = EnergyData(100, 'd')
energy100e = EnergyData(100, 'e')

energy250a = EnergyData(250, 'a')
energy250b = EnergyData(250, 'b')
energy250c = EnergyData(250, 'c')
energy250d = EnergyData(250, 'd')
energy250e = EnergyData(250, 'e')

energy500a = EnergyData(500, 'a')
energy500b = EnergyData(500, 'b')
energy500c = EnergyData(500, 'c')
energy500d = EnergyData(500, 'd')
energy500e = EnergyData(500, 'e')

energy1000a = EnergyData(1000, 'a')
energy1000b = EnergyData(1000, 'b')
energy1000c = EnergyData(1000, 'c')
energy1000d = EnergyData(1000, 'd')
energy1000e = EnergyData(1000, 'e')

energy2000a = EnergyData(2000, 'a')
energy2000b = EnergyData(2000, 'b')
energy2000c = EnergyData(2000, 'c')
energy2000d = EnergyData(2000, 'd')
energy2000e = EnergyData(2000, 'e')

energy5000a = EnergyData(5000, 'a')
energy5000b = EnergyData(5000, 'b')
energy5000c = EnergyData(5000, 'c')
energy5000d = EnergyData(5000, 'd')
energy5000e = EnergyData(5000, 'e')

energy10000a = EnergyData(10000, 'a')
energy10000b = EnergyData(10000, 'b')
energy10000c = EnergyData(10000, 'c')
energy10000d = EnergyData(10000, 'd')
energy10000e = EnergyData(10000, 'e')

timeData = energy100a.time
Nlist = [100, 250, 500, 1000, 2000, 5000, 10000]
Nlist = np.asarray(Nlist)

Elost = []
Elost.append(energy100c.Etot[0] - energy100c.Etot[-1])
Elost.append(energy250c.Etot[0] - energy250c.Etot[-1])
Elost.append(energy500c.Etot[0] - energy500c.Etot[-1])
Elost.append(energy1000c.Etot[0] - energy1000c.Etot[-1])
Elost.append(energy2000c.Etot[0] - energy2000c.Etot[-1])
Elost.append(energy5000c.Etot[0] - energy5000c.Etot[-1])
Elost.append(energy10000c.Etot[0] - energy10000c.Etot[-1])
Elost = np.asarray(Elost)

unbound100 = (energy100a.unbound + energy100b.unbound + energy100c.unbound + energy100d.unbound + energy100e.unbound) / (5*energy100a.N)
unbound250 = (energy250a.unbound + energy250b.unbound + energy250c.unbound + energy250d.unbound + energy250e.unbound) / (5*energy250a.N)
unbound500 = (energy500a.unbound + energy500b.unbound + energy500c.unbound + energy500d.unbound + energy500e.unbound) / (5*energy500a.N)
unbound1000 = (energy1000a.unbound + energy1000b.unbound + energy1000c.unbound + energy1000d.unbound + energy1000e.unbound) / (5*energy1000a.N)
unbound2000 = (energy2000a.unbound + energy2000b.unbound + energy2000c.unbound + energy2000d.unbound + energy2000e.unbound) / (5*energy2000a.N)
unbound5000 = (energy500a.unbound + energy5000b.unbound + energy5000c.unbound + energy5000d.unbound + energy5000e.unbound) / (5*energy5000a.N)
unbound10000 = (energy10000a.unbound + energy10000b.unbound + energy10000c.unbound + energy10000d.unbound + energy10000e.unbound) / (5*energy10000a.N)

plt.figure()
plt.plot(timeData, energy100c.Etot)
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.savefig('plots/plot.energy.total.100.png')

plt.figure()
plt.plot(timeData, energy500c.Etot)
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.savefig('plots/plot.energy.total.500.png')

plt.figure()
plt.plot(timeData, energy1000c.Etot)
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.savefig('plots/plot.energy.total.1000.png')

plt.figure()
plt.plot(timeData, energy10000c.Etot)
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.savefig('plots/plot.energy.total.10000.png')

plt.figure()
plt.plot(timeData, unbound100)
plt.plot(timeData, unbound250)
plt.plot(timeData, unbound500)
plt.plot(timeData, unbound1000)
plt.plot(timeData, unbound2000)
plt.plot(timeData, unbound5000)
plt.plot(timeData, unbound10000)
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('fraction of particles with positive energy')
plt.legend(['N = 100','N = 250','N = 500','N = 1000','N = 2000','N = 5000','N = 10000'], loc=4)
plt.savefig('plots/plot.fraction.positive.energy.png')

plt.figure()
plt.plot(Nlist, Elost)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Energy lost')
plt.savefig('plots/plot.energy.lost.png')

plt.figure()
plt.plot(Nlist, Elost/Nlist)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Energy lost / N')
plt.savefig('plots/plot.energy.lost.N.png')

plt.figure()
plt.plot(timeData, energy1000c.Ep, timeData, energy1000c.Ek, timeData, (2*energy1000c.Ek + energy1000c.Ep))
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.legend(['$E_P$','$E_K$','$2E_K + E_P$'], loc='best')
plt.savefig('plots/plot.energy.virial.1000.png')

plt.figure()
plt.plot(timeData, energy10000c.Ep, timeData, energy10000c.Ek, timeData, (2*energy10000c.Ek + energy10000c.Ep))
plt.xlabel('time [$t_{ff}$]')
plt.ylabel('Energy')
plt.legend(['$E_P$','$E_K$','$2E_K + E_P$'], loc='best')
plt.savefig('plots/plot.energy.virial.10000.png')
