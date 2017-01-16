import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 


filenameRadialBH = 'data/radialBH10000.dat'

try:  
    fileRadialBH = open(filenameRadialBH, 'r') 
except: 
    print '%s not found.' %(filenameRadialBH)
    sys.exit()

radialDistances = []

for line in fileRadialBH:
	radialDistances.append(float(line))

radialDistances = np.asarray(radialDistances)

print 'Average: ', np.average(radialDistances)
print 'Standard deviation: ', np.std(radialDistances)
print max(radialDistances)

m = 500

r = np.linspace(0, 10, m)
n = np.zeros(m)
NFW = np.zeros(100)
n0 = 50.0
r0 = 3.0

rVolume = np.linspace(0, 10, m)
dr = 10./m

for i in range(len(r)-1):
	for j in range(len(radialDistances)):
		if (radialDistances[j] < r[i+1]) and (radialDistances[j] > r[i]):
			n[i+1] +=1

for i in range(len(rVolume)-1):
	rVolume[i+1] = 4/3.0 * np.pi * r[i+1]**3 - rVolume[i]

plt.figure(0)
plt.plot(r[1:]+0.5*dr, n[1:]/(rVolume[1:]))
plt.xscale('log')
plt.yscale('log')
plt.show()
#sys.exit()
sys.exit()
for i in range(100):
	temp = r[i] / r0
	n[i] = n0 / (1.0 + temp**4)
	NFW[i] = n0 / (1.0 + temp * (1.0 + temp*temp))

#plt.plot(r, n, '-r', r, NFW, '-g')
#plt.show()


plt.figure(1)
plt.hist(radialDistances, bins=30, range=[0,10])
plt.xlabel('Distance [ly]')
plt.ylabel('Number')
plt.show()

