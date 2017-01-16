# FYS4150 project 2
# program to print and plot first eigenvalue of omega

import os
import sys
import numpy as np
from matplotlib import pyplot as plt


try:
	n = float(sys.argv[1])
	rhoFinal = float(sys.argv[2])
except:
	print 'Wrong command line arguments.'
	print 'Run as:program n rhoFinal'
	sys.exit()

numberOfElectrons = 2
omega = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0]

datafilename = []

for i in range(len(omega)):
	datafilename.append('data/data.e%d.n%d.rho%g.w%g.dat' %(numberOfElectrons, n, rhoFinal, omega[i]))


firstEigenvalue = []
for i in range(len(omega)):
	# open datafile
	try:
		datafile = open(datafilename[i], 'r')
	except:
		print '%s does not exist.' %(datafilename[i])
		sys.exit()

	linenumber = 0
	for line in datafile:
		if (linenumber == 10):
			words = line.split()
			firstEigenvalue.append(float(words[0]))
			break
		linenumber += 1
	datafile.close()

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


# printing the first eigenvalue for different omega
print 'rhoMax = %f' %(rhoFinal)
for i in range(len(omega)):
	print 'omega: %4g, EV1: %g' %(omega[i], firstEigenvalue[i])

plt.figure(0)
plt.plot(omega, firstEigenvalue, 'b-')
plt.xlabel('Frequenzy $\omega_r$')
plt.ylabel('First eigenvalue $\lambda_0$')

plt.show()