# FYS4150 project 2
# program to compare our code to armadillo

import os
import sys
import numpy as np
from matplotlib import pyplot as plt


# convert command line arguments
try:
	nMin = int(sys.argv[1])
	nMax = int(sys.argv[2])
	rhoFinal = float(sys.argv[3])
except:
	print 'Wrong command line arguments.'
	print 'Run as: program nMin nMax rhoFinal'
	sys.exit()

datafilename = 'data/data.n%d.%d.rhoFinal%g.dat' %(nMin, nMax, rhoFinal)
datafilenameArma = 'data/data.n%d.%d.rhoFinal%g.arma.dat' %(nMin, nMax, rhoFinal)

# open datafile
try:
	datafile = open(datafilename, 'r')
	datafileArma = open(datafilenameArma, 'r')
except:
	print '%s or %s does not exist.' %(datafilename, datafilenameArma)
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
nArma = []
calculationTime = []
calculationTimeArma = []


linenumber = 0

for line in datafile:
	words = line.split()
	if (linenumber > 7):
		n.append(int(words[0]))
		calculationTime.append(float(words[2]))
	linenumber += 1


linenumber = 0

for line in datafileArma:
	words = line.split()
	if (linenumber > 6):
		nArma.append(int(words[0]))
		calculationTimeArma.append(float(words[1]))
	linenumber += 1

datafile.close()
datafileArma.close()

n = np.asarray(n)
nArma = np.asarray(nArma)
calculationTime = np.asarray(calculationTime)
calculationTimeArma = np.asarray(calculationTimeArma)

plt.figure(0)
plt.plot(nArma, calculationTimeArma, 'b-')
plt.xlabel('n (size of matrix)')
plt.ylabel('Calculation time [s]')
plt.savefig('plot.arma.calculationtime.png')


plt.figure(1)
plt.plot(n, calculationTime, 'r-', nArma, calculationTimeArma, 'b-')
plt.xlabel('n (size of matrix)')
plt.ylabel('Calculation time [s]')
plt.xscale('log')
plt.yscale('log')
plt.legend(['Our code', 'Armadillo'], loc='best')
plt.savefig('plot.calculationtime.vs.arma.png')


plt.show()