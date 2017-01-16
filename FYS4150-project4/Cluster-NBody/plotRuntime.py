import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 

N = [100, 250, 500, 1000, 2000, 5000, 10000, 50000, 100000]
BH1 = [2.71, 11.42, 34.5, 109.7, 309.7, 1240, 42.63*50, 203.5*100, 100.3*500]
BH2 = [1.88, 6.85, 19.5, 60.3, 165.0, 692, 22.64*50, 105.2*100, 53.29*500]
BH3 = [2.1, 5.82, 15.6, 45.8, 123.3, 545, 16.12*50, 74.15*100, 38.47*500]
BH4 = [3.7, 6.93, 17.1, 44.8, 109.4, 464, 14.49*50, 62.14*100, 32.68*500]
N1 = [2.72, 15.42, 60.5, 222.4, 23.0*50, 143.7*50, 314.9*100, 1250.0*1000, 2902.0*2500]
N2 = [2.57, 10.12, 37.2, 149.8, 14.8*50, 91.7*50, 189.9*100, 837*1000, 1906*2500]
N3 = [2.74, 9.12, 30.0, 111, 11.6*50, 66.9*50, 139.6*100, 635*1000, 1601*2500]
N4 = [3.17, 8.61, 27.1, 102, 9.57*50, 55.3*50, 113.2*100, 498*1000, 1186*2500]

iterations = 5000.0
Nshort = [100, 250, 500, 1000]
NV1 = np.asarray([2.72, 15.42, 60.5, 222.4])
NRK1 = np.asarray([6.98, 37.09, 151.6, 614.7])


plt.figure(0)
plt.plot(N, BH1, N, BH2, N, BH3, N, BH4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Runtime [s]')
plt.legend(['1 thread','2 threads','3 threads','4 threads'], loc='best')
plt.savefig('plots/plot.runtime.BH.png')

plt.figure(1)
plt.plot(N, N1, N, N2, N, N3, N, N4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Runtime [s]')
plt.legend(['1 thread','2 threads','3 threads','4 threads'], loc='best')
plt.savefig('plots/plot.runtime.NN.png')


plt.figure(2)
plt.plot(N, N1, N, BH1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Runtime [s]')
plt.legend(['Direct','Barnes-Hut'], loc='best')
plt.savefig('plots/plot.runtime.comparison.png')


plt.figure(3)
plt.plot(Nshort, NV1/iterations, Nshort, NRK1/iterations)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Runtime per iteration [s]')
plt.legend(['Velocity Verlet','Runge-Kutta 4'], loc='best')
plt.savefig('plots/plot.runtime.iteration.png')


#plt.show()
