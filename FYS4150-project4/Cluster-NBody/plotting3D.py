# FYS415 Cluster
# plotting and animating data

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

filename = 'data/output.dat'

try:
    fileposition = open(filename, 'r')
except:
    print 'failed to open %s' %(filename)
N = 100
sideLength = 40
steps = int(2 / 0.001 + 1)


x = np.zeros((N, steps))
y = np.zeros((N, steps))
z = np.zeros((N, steps))

counter = 0
for line in fileposition:
    words = line.split()
    for i in range(N):
        x[i, counter] = float(words[i])
        y[i, counter] = float(words[i+1])
        z[i, counter] = float(words[i+2])
    counter += 1

'''
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
for i in range(100):
    plt.plot(x[i,:], y[i,:], z[i,:])
ax.set_xlim3d([-10, 10])
ax.set_ylim3d([-10, 10])
ax.set_zlim3d([-10, 10])
#plt.ylim([-30, 30])
plt.show()
'''
class AnimatedScatter3D(object):
    def __init__(self, points, x, y, z, steps, sideLength):
        self.N = points
        self.x = x
        self.y = y
        self.z = z
        self.steps = steps
        self.sideLength = sideLength
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.anim = animation.FuncAnimation(self.fig, self.animateUpdate, interval=100, init_func=self.setupPlot, blit=True)
    
    def animateUpdate(self, i):
        update = []
        #print type((x[0,i], y[0,i], z[0,i]))
        for j in range(N):
            update.append([x[j,i], y[j,i], z[j,i]])
        #update = np.asarray(update)
        print type(update)
        #sys.exit()
        self.scat._offsets3d(update)

        #self.scat._offsets3d(np.ma.ravel(x[:,0]), np.ma.ravel(y[:,0]), np.ma.ravel(z[:,0]))
        plt.draw()
        return self.scat,

    def setupPlot(self):
        self.scat = self.ax.scatter(x[:,0], y[:,0], z[:,0], animated=True)
        self.ax.set_xlim3d(-sideLength/2.0, sideLength/2.0)
        self.ax.set_ylim3d(-sideLength/2.0, sideLength/2.0)
        self.ax.set_zlim3d(-sideLength/2.0, sideLength/2.0)
        return self.scat,

    def show(self):
        plt.show()

a = AnimatedScatter3D(N, x, y, z, steps, sideLength)
a.show()
