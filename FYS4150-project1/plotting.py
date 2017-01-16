# FYS4150-project1
# program to plot data

import os
import sys
import numpy as np
from matplotlib import pyplot as plt


# convert command line argument to int
try:
	n = int(sys.argv[1])
except:
	print 'Wrong command line arguments.'
	print 'Run as: program integer'
	sys.exit()

datafilename = 'data.n%d.dat' %(n)

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

h = 0

linenumber = 0
# first line contains n and h
for line in datafile:
	words = line.split()
	#n = int(words[0])
	h = float(words[1])
	linenumber += 1
	break

# for storage
v = np.zeros(n+2)
x = np.linspace(0, 1, n+2)
relError = np.zeros(n+2)
vLU = np.zeros(n+2)
continsvLU = False		# boolean to keep track of wether LU data is saved to file

# reading values from file and adding to arrays
for line in datafile:
	if (linenumber > 0):
		words = line.split()
		v[linenumber-1] = float(words[0])
		relError[linenumber-1] = float(words[1])
		# contains vLU
		if ((len(words) > 2) or (continsvLU)):
			vLU[linenumber-1] = float(words[2])
			continsvLU = True
		linenumber += 1

# close the file
datafile.close()

# finding values for the exact solution
values = 10001
u_closed_form = np.zeros(values)
x_values = np.linspace(0, 1, values)
for i in range(values):
	u_closed_form[i] = 1 - (1 - np.exp(-10)) * x_values[i] - np.exp(-10 * x_values[i])

# plotting numerical vs closed form solution
plt.figure(0)
plt.plot(x_values, u_closed_form, x, v)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.legend(['closed form solution', 'numerical $n = 10^{%g}$' %(np.log10(n))], loc=0)
plt.savefig('plot.n%d.vs.png' %(n))

# plotting the relative error
plt.figure(1)
plt.plot(x, relError)
plt.xlabel('$x$')
plt.ylabel('relError')
plt.savefig('plot.n%d.error.png' %(n))

# plotting v vs vLU
if (continsvLU):
	plt.figure(2)
	plt.plot(x, v, '-g', x, vLU, '-r')
	plt.xlabel('$x$')
	plt.ylabel('$u(x)$')
	plt.legend(['numerical $n = 10^{%g}$' %(np.log10(n)), 'LU decomposition'], loc=0)
	plt.savefig('plot.n%d.vsLU.png' %(n))

print 'plots saved to %s' %(plotpath)
plt.show()
