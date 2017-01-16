import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 


filenameEnergyBH = 'data/compareenergyBH1000.dat'
filenameEnergyNtoN = 'data/compareenergyNtoN1000.dat'

try:  
    fileEnergyBH = open(filenameEnergyBH, 'r') 
except: 
    print '%s not found.' %(filenameEnergyBH)
    sys.exit()

try:  
    fileEnergyNtoN = open(filenameEnergyNtoN, 'r') 
except: 
    print '%s not found.' %(filenameEnergyNyoN)
    sys.exit() 

timeBH = []
energyKineticBH = []
energyPotentialBH = []
energyTotalBH = []

for line in fileEnergyBH:
    words = line.split()
    timeBH.append(float(words[0]))
    energyKineticBH.append(float(words[1]))
    energyPotentialBH.append(float(words[2]))

timeBH = np.asarray(timeBH)
energyKineticBH = np.asarray(energyKineticBH)
energyPotentialBH = np.asarray(energyPotentialBH)
energyTotalBH = energyKineticBH + energyPotentialBH

timeNtoN = []
energyKineticNtoN = []
energyPotentialNtoN = []
energyTotalNtoN = []

for line in fileEnergyNtoN:
    words = line.split()
    timeNtoN.append(float(words[0]))
    energyKineticNtoN.append(float(words[1]))
    energyPotentialNtoN.append(float(words[2]))

timeNtoN = np.asarray(timeNtoN)
energyKineticNtoN = np.asarray(energyKineticNtoN)
energyPotentialNtoN = np.asarray(energyPotentialNtoN)
energyTotalNtoN = energyKineticNtoN + energyPotentialNtoN


plt.figure(0)
plt.plot(timeBH, energyKineticBH, timeBH, energyPotentialBH, timeBH, energyTotalBH)
#plt.plot(timeBH, energyPotentialBH)
plt.legend(['$E_k$','$E_p$','$E_{tot}$'])
#plt.title('BarnesHut')
plt.xlabel('Time [$\\tau_{crunch}$]')
plt.ylabel('Energy')
plt.savefig('plots/compareBH.png')

plt.figure(1)
plt.plot(timeNtoN, energyKineticNtoN, timeNtoN, energyPotentialNtoN, timeNtoN, energyTotalNtoN)
#plt.plot(timeNtoN, energyPotentialNtoN)
plt.legend(['$E_k$','$E_p$','$E_{tot}$'])
#plt.title('NtoN')
plt.xlabel('Time [$\\tau_{crunch}$]')
plt.ylabel('Energy')
plt.savefig('plots/compareNtoN.png')

plt.figure(2)
plt.plot(timeNtoN, energyTotalNtoN, timeBH, energyTotalBH)
plt.legend(['Direct','Barnes-Hut'])
plt.xlabel('Time [$\\tau_{crunch}$]')
plt.ylabel('Energy')
plt.savefig('plots/compareEnergy.png')

#plt.figure(2)
#plt.plot(timeBH, energyPotentialBH, timeNtoN, energyPotentialNtoN)
#plt.legend(['BH','NtoN'])
#plt.xlabel('Time [$\\tau_{crunch}$]')
#plt.ylabel('Energy')
#plt.savefig('plots/energyCompare200.png')

#plt.show()
