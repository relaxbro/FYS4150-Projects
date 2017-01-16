import os 
import sys 
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import rcParams 
rcParams.update({'figure.autolayout': True}) 

def gauss(x, mu, sigma):
    return 1.0 / (sigma * np.sqrt(2*np.pi)) * np.exp(- (x-mu)**2 / (2*sigma**2))

try: 
    N = int(sys.argv[1])
    letter = sys.argv[2]
except:
    print 'Missing command line arguments. program N letter'
    sys.exit()


#N = 10000
#letter = 'c'
filenameRadialBH = 'data/radialBH%d%s.dat' %(N, letter)
R0 = 20

try:  
    fileRadialBH = open(filenameRadialBH, 'r') 
except: 
    print '%s not found.' %(filenameRadialBH)
    sys.exit()

CoM = np.zeros(3)
x = []
y = []
z = []

linenumber = 0
for line in fileRadialBH:
 	if (linenumber == 0):
 		words = line.split()
 		CoM[0] = float(words[0])
 		CoM[1] = float(words[1])
 		CoM[2] = float(words[2])
 		linenumber += 1
 		break

for line in fileRadialBH:
	words = line.split()
 	x.append(float(words[0]))
 	y.append(float(words[1]))
 	z.append(float(words[2]))

x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)

xmu = np.average(x)
xsigma = np.std(x)
ymu = np.average(y)
ysigma = np.std(y)
zmu = np.average(z)
zsigma = np.std(z)

print 'Center of mass: ', CoM

xbin = np.round(CoM[0])
ybin = np.round(CoM[1])
zbin = np.round(CoM[2])

xdata = []
ydata = []
zdata = []

for i in range(len(x)):
    if (abs(x[i] - xbin) < 4):
        if (abs(y[i] - ybin) < 4):
            if (abs(z[i] - zbin) < 4):
                xdata.append(x[i])
                ydata.append(y[i])
                zdata.append(z[i])

xmu = np.average(xdata)
xsigma = 0.75*np.std(xdata)
ymu = np.average(ydata)
ysigma = 0.75*np.std(ydata)
zmu = np.average(zdata)
zsigma = 0.75*np.std(zdata)

'''
plt.figure(1)
#plt.hist(x, bins=100, range=[-20,20])
count, xbins, ignored = plt.hist(xdata, bins=50, range=[xbin-4, xbin+4])
xga = gauss(xbins, xmu, xsigma)
plt.plot(xbins, max(count) * xga / max(xga), linewidth=2, color='r')
plt.xlabel('Distance [ly]')
plt.ylabel('Number per bin')
#plt.savefig('plots/plot.radial.dist.hist.10000x.png')

plt.figure(2)
#plt.hist(y, bins=100, range=[-20,20])
count, ybins, ignored = plt.hist(ydata, bins=50, range=[ybin-4, ybin+4])
yga = gauss(ybins, ymu, ysigma)
plt.plot(ybins, max(count) * yga / max(yga), linewidth=2, color='r')
plt.xlabel('Distance [ly]')
plt.ylabel('Number per bin')
#plt.savefig('plots/plot.radial.dist.hist.10000y.png')

plt.figure(3)
#plt.hist(z, bins=100, range=[-20,20])
count, zbins, ignored = plt.hist(zdata, bins=50, range=[zbin-4, zbin+4])
zga = gauss(zbins, zmu, zsigma)
plt.plot(zbins, max(count) * zga / max(zga), linewidth=2, color='r')
plt.xlabel('Distance [ly]')
plt.ylabel('Number per bin')
#plt.savefig('plots/plot.radial.dist.hist.10000z.png')
'''



RtoCoM = np.zeros(len(x))

for i in range(len(x)):
    RtoCoM[i] = np.sqrt( (x[i] - CoM[0])**2 + (y[i] - CoM[1])**2 + (z[i] - CoM[2])**2)

print 'Average: ', np.average(RtoCoM)
print 'Standard deviation: ', np.std(RtoCoM)
print 'Max distance: ', max(RtoCoM)
sys.exit()



m = 500
rmax = 10

r = np.linspace(0, rmax, m)
n = np.zeros(m)
nfit = np.zeros(m)
NFW = np.zeros(m)
V0 = 4.0 / 3 * np.pi * R0**3
n0 = N*N / V0
#n0 = 0.01 * N
r0 = 1.0 / 3 *R0 * N**(-1.0/3)

rVolume = np.linspace(0, rmax, m)
dr = float(rmax)/m

for i in range(len(r)-1):
    for j in range(len(RtoCoM)):
        if (RtoCoM[j] < r[i+1]) and (RtoCoM[j] > r[i]):
            n[i+1] +=1

for i in range(len(rVolume)-1):
    rVolume[i+1] = 4/3.0 * np.pi * r[i+1]**3 - rVolume[i]

for i in range(m):
    temp = r[i] / r0
    nfit[i] = n0 / (1.0 + temp**4)
    NFW[i] = n0 / (1.0 + temp * (1.0 + temp*temp))

plt.figure()
plt.plot(r[1:]+0.5*dr, n[1:]/(rVolume[1:]), 'b')
plt.plot(r[1:], nfit[1:], 'r')
#plt.plot(r[1:], NFW[1:], 'g')
plt.xlim(0.01, 10)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [ly]')
plt.ylabel('n(r) [ly$^{-3}$]')
plt.legend(['N = %d' %(N),'$1/(1 + (r/r_0)^4)$'], loc='best')
#plt.savefig('plots/plot.radial.dist.%d%s.png' %(N, letter))
