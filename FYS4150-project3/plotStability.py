import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


# base filename as commandline argument
try:
    filenameBase = sys.argv[1]
except:
    print 'Run as: program filenameBase'
    sys.exit()


# filenames
filenamePosition = 'data/position.' + filenameBase + '.dat'

try:
    filePosition = open(filenamePosition, 'r')
except:
    print '%s not found.' %(filenamePosition)
    sys.exit()


groundpath = sys.path[0]
plotpath = sys.path[0] + '/plots/'

# change/create plotpath
try:
    os.chdir(plotpath)
except:
    try:
        os.mkdir('plots')
        os.chdir(plotpath)
    except:
        'Failed to create subfolder /plots/. Please create it manually.'
        sys.exit()


timeEnd = 0
dt = 0
bodies = []
mass = []

# read beginning of first file to find initial data
linenumber = 0
for line in filePosition:
    words = line.split()
    if (linenumber == 1):
        timeEnd = float(words[1])
    if (linenumber == 2):
        dt = float(words[1])
    if (linenumber == 3):
        numberOfBodies = int(len(words) / 2)
        for i in range(numberOfBodies):
            bodies.append(words[1 + 2*i])
            mass.append(float(words[2 + 2*i]))
        linenumber += 1
        break
    linenumber += 1

timeSteps = int(timeEnd / dt) + 1


mass = np.asarray(mass)
time = np.zeros(timeSteps)
x = np.zeros((numberOfBodies, timeSteps))
y = np.zeros((numberOfBodies, timeSteps))
z = np.zeros((numberOfBodies, timeSteps))


# read file for positions
index = 0
for line in filePosition:
    if (linenumber > 4):
        words = line.split()
        if (index == timeSteps):
            break
        time[index] = float(words[0])
        for i in range(numberOfBodies):
            x[i, index] = float(words[1 + (3*i)])
            y[i, index] = float(words[2 + (3*i)])
            z[i, index] = float(words[3 + (3*i)])
        index += 1      
    linenumber += 1
filePosition.close()

# findig maximum value for x and y
xmax = 0
ymax = 0
for i in range(numberOfBodies):
    for j in range(timeSteps):
        if (abs(x[i,j]) > xmax):
            xmax = abs(x[i,j])
        if (abs(y[i,j]) > ymax):
            ymax = abs(y[i,j])
if (xmax > ymax):
    ymax = xmax
else:
    xmax = ymax

plt.figure(0)
#ax = plt.axes(xlim=(-1.6*max(x[2,:]), 1.6*max(x[2,:])), ylim=(-1.2*max(y[2,:]), 1.2*max(y[2,:])))
plt.plot([0], [0], 'o')
plt.plot(x[1,:], y[1,:])
plt.xlim(-1.6, 1.6)
plt.ylim(-1.2, 1.2)
#plt.axis('equal')
plt.xlabel('Distance [AU]')
plt.ylabel('Distance [AU]')
plt.legend(bodies)
savefilepos = "integration." + filenameBase + '.pos.png'
plt.savefig(savefilepos)
#plt.show()