# FYS4150-project1
# program to plot number of iterations

import os
import sys
import numpy as np
from matplotlib import pyplot as plt


# convert command line argument
try:
	nMin = int(sys.argv[1])
	nMax = int(sys.argv[2])
	rhoFinal = float(sys.argv[3])
except:
	print 'Wrong command line arguments.'
	print 'Run as: program nMin nMax rhoFinal'
	sys.exit()

datafilename = 'data/data.n%d.%d.rhoFinal%g.dat' %(nMin, nMax, rhoFinal)

# open datafile
try:
	datafile = open(datafilename, 'r')
except:
	print '%s does not exist.' %(datafilename)
	sys.exit()

plotpath = sys.path[0] + '/plots/'

# changing working directory for saving plots
try:
	os.chdir(plotpath)
except:
	try:
		os.mkdir('plots')
		os.chdir(plotpath)
	except:
		print 'Failed to create subfolder /plots/. Please create it manually.'


n = []
iterations = []
calculationTime = []
tolerance = 0


linenumber = 0

for line in datafile:
	words = line.split()
	if (linenumber == 5):
		tolerance = words[1]
	elif (linenumber > 7):
		n.append(int(words[0]))
		iterations.append(int(words[1]))
		calculationTime.append(float(words[2]))
	linenumber += 1

# close the file
datafile.close()

n = np.asarray(n)
iterations = np.asarray(iterations)
calculationTime = np.asarray(calculationTime)

# finding polynomial fit
polyarr = np.polyfit(n, iterations, 2)
p = np.poly1d(polyarr)

print p


narr = np.linspace(0, nMax+100, 101)

plt.figure(0)
plt.plot(n, calculationTime, 'r-')
plt.xlabel('n (size of matrix)')
plt.ylabel('Calculation time [s]')
plt.savefig('plot.calculationtime.png')


plt.figure(1)
plt.plot(n, iterations, 'b-')
plt.xlabel('n (size of matrix)')
plt.ylabel('iterations')
plt.savefig('plot.iterations.png')


plt.figure(2)
plt.plot(n, iterations, 'b-', narr, p(narr), 'r-')
plt.xlabel('n (size of matrix)')
plt.ylabel('iterations')
plt.legend(['data', 'best fit'], loc='best')
plt.savefig('plto.iterations.bestfit.png')

plt.show()
