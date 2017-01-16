# FYS4150 project 3
# program to convert and write initial data to files

import os
import sys
import numpy as np

# inputfile in SI units
try:   
    datafile = open('initialdata/solarsystemSI.txt', 'r')
except:
    print 'missing file'
    sys.exit()

# empty lists
numberOfBodies = 0
bodies = []
mass = []
x = []
y = []
z = []
vx = []
vy = []
vz = []

# read all data from file
linenumber = 0
for line in datafile:
    if (linenumber > 5):
        words = line.split()
        bodies.append(words[0])
        x.append(float(words[1]))
        y.append(float(words[2]))
        z.append(float(words[3]))
        vx.append(float(words[4]))
        vy.append(float(words[5]))
        vz.append(float(words[6]))
        mass.append(float(words[7]))
        numberOfBodies += 1
    linenumber += 1
datafile.close()

# convert to arrays
x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)
vx = np.asarray(vx)
vy = np.asarray(vy)
vz = np.asarray(vz)
mass = np.asarray(mass)

# conversion factors
mSun = 1.99e30
meterToAU = 1.496e11
metersPerSecontToAUPerYear = 4741

# convert to AU yr m_sun units
mass = mass / mSun
x = x / (meterToAU * 2 * np.pi)
y = y / (meterToAU * 2 * np.pi)
z = z / (meterToAU * 2 * np.pi)
vx = vx / metersPerSecontToAUPerYear
vy = vy / metersPerSecontToAUPerYear
vz = vz / metersPerSecontToAUPerYear


# find suns position for center of mass to be inn 0,0
mtimesR = 0
for i in range(1, numberOfBodies, 1):
    mtimesR += x[i] * mass[i]
x[0] = - mtimesR / mass[0]

# find suns momentum for total momentum to be zero
momentum = 0
for i in range(1, numberOfBodies, 1):
    momentum += mass[i] * vy[i]
vy[0] = -momentum / mass[0]

# write whole solar system to file
outfile = open('initialdata/solarsystemAU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
outfile.close()

# write file for sun + earth, sun in placed in 0,0
outfile = open('initialdata/sunearthAU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
outfile.close()

# write file for sun + earth, sun in placed in 0,0 and earth escaping
outfile = open('initialdata/sunearthescapeAU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], 8.9, vz[i], mass[i]))
outfile.close()

# write file for sun + earth + moon, sun in 0,0
outfile = open('initialdata/sunearthmoonAU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth' or bodies[i] == 'moon'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
outfile.close()

# write file for sun + earth + jupiter with sun at 0,0
outfile = open('initialdata/sunearthjupiterAU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth' or bodies[i] == 'jupiter'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
outfile.close()

# write file with sun + earth + jupiter with sun at 0,0 and jupiter mass increased by 10
outfile = open('initialdata/sunearthjupiter10AU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
    elif (bodies[i] == 'jupiter'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], 10*mass[i]))
outfile.close()

# write file with sun + earth + jupiter with sun at 0,0 and jupiter mass increased by 1000
outfile = open('initialdata/sunearthjupiter1000AU.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')
for i in range(numberOfBodies):
    if (bodies[i] == 'sun'):
        outfile.write('%s 0 0 0 0 0 0 %e\n' %(bodies[i], mass[i]))
    elif (bodies[i] == 'earth'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
    elif (bodies[i] == 'jupiter'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], 1000*mass[i]))
outfile.close()

# write file for sun + earth + jupiter with total momentum 0
outfile = open('initialdata/sunearthjupiterAUmotion.txt', 'w')
outfile.write('units\ndistance: AU\nvelocity: AU/yr\nmass: solarmasses\n')
outfile.write('data\n')
outfile.write('name -- x -- y -- z -- vx -- vy -- vz -- mass\n')

mtimesR = 0
for i in range(1, numberOfBodies, 1):
    if (bodies[i] == 'earth' or bodies[i] == 'jupiter'):
        mtimesR += x[i] * mass[i]
x[0] = - mtimesR / mass[0]

momentum = 0
for i in range(1, numberOfBodies, 1):
    if (bodies[i] == 'earth' or bodies[i] == 'jupiter'):
        momentum += mass[i] * vy[i]
vy[0] = -momentum / mass[0]

for i in range(numberOfBodies):
    if (bodies[i] == 'earth' or bodies[i] == 'jupiter' or bodies[i] == 'sun'):
        outfile.write('%s %e %e %e %e %e %e %e\n' %(bodies[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], mass[i]))
outfile.close()
