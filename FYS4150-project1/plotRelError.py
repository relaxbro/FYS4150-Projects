# FYS4150-project1
# program to plot relative error

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

n = []
# convert command line argument to int
if (len(sys.argv) > 1):
	try:
		for i in range(1, len(sys.argv)):
			n.append(int(sys.argv[i]))
		n = np.asarray(n)	
	except:
		print 'Wrong command line arguments.'
		print 'Run as: program integer integer ...'
		n = np.asarray([10, 100, 1000, 10000, 100000, 1000000])
else:
	print 'Missing command line arguments, trying 10 100 1000 10000 100000 1000000.'
	n = np.asarray([10, 100, 1000, 10000, 100000, 1000000])

datafilenames = []

# creating list of filenames
for i in n:
	datafilenames.append('data.n%d.dat' %(i))

datafiles = []

# open all datafiles
for datafilename in datafilenames:
	try:
		datafiles.append(open(datafilename, 'r'))
	except:
		print '%s does not exist.' %(datafilename)
		print 'chose existing data sets'
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

# for storing
h = np.zeros(len(n))
relErrorMax = np.zeros(len(n))
v = []
relError = []
x = []

for i in n:
	v.append(np.zeros(i+2))
	relError.append(np.zeros(i+2))
	x.append(np.linspace(0, 1, i+2))

# convert to arrays
v = np.asarray(v)
relError = np.asarray(relError)
x = np.asarray(x)

for i in range(len(datafiles)):
	linenumber = 0
	# first line contains n, h and max relative error
	for line in datafiles[i]:
		words = line.split()
		h[i] = float(words[1])
		relErrorMax[i] = float(words[2])
		linenumber += 1
		break

	for line in datafiles[i]:
		if (linenumber > 0):
			words = line.split()
			v[i][linenumber-1] = float(words[0])
			relError[i][linenumber-1] = float(words[1])
			linenumber += 1

# close all files
for datafile in datafiles:
	datafile.close()


# plotting the relative error
fig = plt.figure(0)
ax = plt.subplot(111)
for i in range(len(x)):
	line, = ax.plot(x[i], relError[i], label='$n = 10^{%g}$' %(np.log10(n[i])))
plt.xlabel('$x$')
plt.ylabel('relError')
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=5)
plt.savefig('plot.all.error.png')

# plotting the maximum relative error as a function of log10(h)
plt.figure(1)
plt.plot(h, relErrorMax)
plt.xscale('log')
plt.xlabel('$log_{10}(h)$')
plt.ylabel('max relError')
plt.legend(['max $\epsilon$'], loc=2)
plt.savefig('plot.all.errormax.png')

print 'plots saved to %s' %(plotpath)
plt.show()