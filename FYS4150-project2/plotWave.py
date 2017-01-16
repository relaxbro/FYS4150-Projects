# FYS4150 project 2
# program to plot wavefunction

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

# values to find out datafiles
n = 200
numberOfElectrons = 2
rhoFinal = 5.0
rho = np.linspace(0, rhoFinal, n)
#omega = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0]
omega = [0.01, 0.5, 1.0, 5.0]

#create file names
datafilename = []
for i in range(len(omega)):
	datafilename.append('data/data.e%d.n%d.rho%g.w%g.dat' %(numberOfElectrons, n, rhoFinal, omega[i]))

# array to contain all eigenvectors
eigenvector1 = np.zeros((len(omega), n))
eigenvector2 = np.zeros((len(omega), n))
eigenvector3 = np.zeros((len(omega), n))

for i in range(len(omega)):
	# open datafile
	try:
		datafile = open(datafilename[i], 'r')
	except:
		print '%s does not exist.' %(datafilename[i])
		sys.exit()

	linenumber = 0
	counter = 0
	for line in datafile:
		if (linenumber > 9):
			words = line.split()
			eigenvector1[i, counter] = float(words[1])
			eigenvector2[i, counter] = float(words[2])
			eigenvector3[i, counter] = float(words[3])
			counter += 1
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


plt.figure(0)
# array containing legend names
legendarr = []
# plot all eigenvectors
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector1[i,:], eigenvector1[i,:])
	eigenvector1[i,:] =  eigenvector1[i,:] / np.sqrt(norm)
	plt.plot(rho, eigenvector1[i,:])
	legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Wavefuntion $\psi$')
plt.title('Ground state')
plt.legend(legendarr)
plt.savefig('plot.wave1.png')

plt.figure(1)
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector2[i,:], eigenvector2[i,:])
	eigenvector2[i,:] =  eigenvector2[i,:] / np.sqrt(norm)
	plt.plot(rho, eigenvector2[i,:])
	#legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Wavefuntion $\psi$')
plt.title('First excited state')
plt.legend(legendarr)
plt.savefig('plot.wave2.png')

plt.figure(2)
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector3[i,:], eigenvector3[i,:])
	eigenvector3[i,:] =  eigenvector3[i,:] / np.sqrt(norm)
	plt.plot(rho, eigenvector3[i,:])
	#legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Wavefuntion $\psi$')
plt.title('Second excited state')
plt.legend(legendarr, loc=4)
plt.savefig('plot.wave3.png')

plt.figure(3)
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector1[i,:], eigenvector1[i,:])
	eigenvector1[i,:] =  eigenvector1[i,:] / np.sqrt(norm)
	plt.plot(rho, np.sqrt(eigenvector1[i,:]**2))
	#legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Probability $|\psi|$')
plt.title('Ground state')
plt.legend(legendarr)
plt.savefig('plot.wave1p.png')

plt.figure(4)
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector2[i,:], eigenvector2[i,:])
	eigenvector2[i,:] =  eigenvector2[i,:] / np.sqrt(norm)
	plt.plot(rho, np.sqrt(eigenvector2[i,:]**2))
	#legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Probability $|\psi|$')
plt.title('First excited state')
plt.legend(legendarr)
plt.savefig('plot.wave2p.png')

plt.figure(5)
for i in range(len(omega)):
	# normalization factor
	norm = np.inner(eigenvector3[i,:], eigenvector3[i,:])
	eigenvector3[i,:] =  eigenvector3[i,:] / np.sqrt(norm)
	plt.plot(rho, np.sqrt(eigenvector3[i,:]**2))
	#legendarr.append('$\omega_r = %g$' %(omega[i]))
plt.xlabel('$\\rho$')
plt.ylabel('Probability $|\psi|$')
plt.title('Second excited state')
plt.legend(legendarr)
plt.savefig('plot.wave3p.png')

plt.show()