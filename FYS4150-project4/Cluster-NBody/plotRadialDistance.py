import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 

N = np.asarray([100, 250, 500, 1000, 2000, 5000, 10000])
Ravg = np.zeros(len(N))
Rstd = np.zeros(len(N))

Ravg[0] = 10.99 + 8.46 + 9.99 + 10.08 + 9.8
Ravg[1] = 9.6 + 12.6 + 7.8 + 6.8 + 7.7
Ravg[2] = 8.1 + 7.7 + 9.5 + 7.6 + 8.3
Ravg[3] = 5.3 + 7.9 + 4.9 + 7.3 + 7.7
Ravg[4] = 4.6 + 4.7 + 6.6 + 4.9 + 5.2
Ravg[5] = 4.7 + 3.9 + 4.7 + 3.1 + 4.2
Ravg[6] = 2.6 + 3.1 + 3.1 + 3.3 + 2.8
Ravg /= 5.0

Rstd[0] = 15.29 + 8.98 + 14.39 + 15.21 + 12.8
Rstd[1] = 15.2 + 16.3 + 11.4 + 10.1 + 11.6
Rstd[2] = 14 + 13.4 + 16.4 + 13.3 + 13.9
Rstd[3] = 10.5 + 13.8 + 10.1 + 13.8 + 12.3
Rstd[4] = 9.8 + 10.7 + 12.7 + 10.2 + 11
Rstd[5] = 11 + 9.3 + 9.7 + 8.1 + 10
Rstd[6] = 7.5 + 8.3 + 8.1 + 8.2 + 7.8
Rstd /= 5.0

plt.figure()
plt.plot(N, Ravg)
plt.ylim(0, 10)
plt.xlabel('N')
plt.ylabel('Average distance to center of mass [ly]')
plt.savefig('plots/plot.radial.distance.average.png')

plt.figure()
plt.plot(N, Rstd)
plt.ylim(0, 15)
plt.xlabel('N')
plt.ylabel('Standard deviation [ly]')
plt.savefig('plots/plot.radial.distance.sigma.png')


