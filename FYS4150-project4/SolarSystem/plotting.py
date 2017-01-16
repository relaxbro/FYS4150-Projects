# FYS415 project 3
# plotting and animating data

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
filenameVelocity = 'data/velocity.' + filenameBase + '.dat'
filenameForce = 'data/force.' + filenameBase + '.dat'
filenameEnergy = 'data/energy.' + filenameBase + '.dat'
filenameAngularMomentum = 'data/angular.momentum.' + filenameBase +'.dat'

# open all files
try:
	filePosition = open(filenamePosition, 'r')
except:
	print '%s not found.' %(filenamePosition)
	sys.exit()

try:
	fileVelocity = open(filenameVelocity, 'r')
except:
	print '%s not found.' %(filenameVelocity)
	sys.exit()

try:
	fileForce = open(filenameForce, 'r')
except:
	print '%s not found.' %(filenameForce)
	sys.exit()

try: 
	fileEnergy = open(filenameEnergy, 'r')
except:
	print '%s not found.' %(filenameEnergy)
	sys.exit()

try:
	fileAngularMomentum = open(filenameAngularMomentum, 'r')
except:
	print '%s not found.' %(filenameAngularMomentum)
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
dx = np.zeros((numberOfBodies, timeSteps))
dy = np.zeros((numberOfBodies, timeSteps))
dz = np.zeros((numberOfBodies, timeSteps))
Fx = np.zeros((numberOfBodies, timeSteps))
Fy = np.zeros((numberOfBodies, timeSteps))
Fz = np.zeros((numberOfBodies, timeSteps))
Ek = np.zeros(timeSteps)
Ep = np.zeros(timeSteps)
Etot = np.zeros(timeSteps)
Lx = np.zeros(timeSteps)
Ly = np.zeros(timeSteps)
Lz = np.zeros(timeSteps)
Ltot = np.zeros(timeSteps)

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

# read file for velocities
linenumber = 0
index = 0
for line in fileVelocity:
	if (linenumber > 4):
		words = line.split()
		if (index == timeSteps):
			break
		for i in range(numberOfBodies):
			dx[i, index] = float(words[1 + (3*i)])
			dy[i, index] = float(words[2 + (3*i)])
			dz[i, index] = float(words[3 + (3*i)])
			#tre = 0
		index += 1	
	linenumber += 1
fileVelocity.close()

# red file for forces
linenumber = 0
index = 0
for line in fileForce:
	if (linenumber > 4):
		words = line.split()
		if (index == timeSteps):
			break
		for i in range(numberOfBodies):
			Fx[i, index] = float(words[1 + (3*i)])
			Fy[i, index] = float(words[2 + (3*1)])
			Fz[i, index] = float(words[3 + (3*i)])
		index += 1
	linenumber += 1
fileForce.close()

# read file for energies
linenumber = 0
index = 0
for line in fileEnergy:
	if (linenumber > 4):
		words = line.split()
		if (index == timeSteps):
			break
		Ek[index] = float(words[1])
		Ep[index] = float(words[2])
		Etot[index] = float(words[3])
		index += 1
	linenumber += 1
fileForce.close()

# read file for angular momentum
linenumber = 0
index = 0
for line in fileAngularMomentum:
	if (linenumber > 4):
		words = line.split()
		if (index == timeSteps):
			break
		Lx[index] = float(words[1])
		Ly[index] = float(words[2])
		Lz[index] = float(words[3])
		index += 1
	linenumber += 1
fileAngularMomentum.close()

# initial function for animate, not used
def init():
	line.set_data([], [])
	return line,

# function for line animation
def animate(i):
	thisx = []
	thisy = []
	for j in range(numberOfBodies):
		thisx.append(x[j,i])
		thisy.append(y[j,i])	
	line.set_data(thisx, thisy)
	return line,

# function for scatter animation
def animatePoints(i, fig, scat, line):
	update = []
	for j in range(numberOfBodies):
		if (i == 0):
			update.append((xmax*0.85, ymax*(1.02-0.1*j)))
		else:
			update.append((x[j,i], y[j,i]))
	line.set_data([-xmax*0.8, xmax*(2*float(i)/timeSteps - 0.8)], [-ymax*1.07, -ymax*1.07])
	scat.set_offsets(update)
	return scat, line,# annotation

# function for scatter animation with framediver
def animatePointsSave(i, frameDivider, fig, scat, line):
	update = []
	i = frameDivider * i
	for j in range(numberOfBodies):
		if (i == 0):
			update.append((xmax*0.85, ymax*(1.02-0.1*j)))
		else:
			update.append((x[j,i], y[j,i]))
	line.set_data([-xmax*0.8, xmax*(2*float(i)/timeSteps - 0.8)], [-ymax*1.07, -ymax*1.07])
	scat.set_offsets(update)
	return scat, line,

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

# finding size of scatter after mass
massSorted = np.sort(mass)
dotSize = np.zeros(numberOfBodies)
for i in range(numberOfBodies):
	for j in range(numberOfBodies):
		if (mass[i] == massSorted[j]):
			if (bodies[i] == 'moon'):
				dotSize[i] = 0.1
			else: 
				dotSize[i] = (0.8) ** (abs(j - numberOfBodies+1))
			break


# animate if wanted
if (len(sys.argv) > 2 and sys.argv[2] == 'animate'):

	# check if save animation and frames wanted
	saveAnimation = False
	if (len(sys.argv) > 3 and sys.argv[3] == 'save'):
		saveAnimation = True
		numberOfFrames = timeSteps
	elif (len(sys.argv) > 3):
		numberOfFrames = int(sys.argv[3])
	else:
		numberOfFrames = timeSteps
	
	# creathe plot
	fig = plt.figure(0)
	ax = plt.axes(xlim=(-1.6*xmax, 1.6*xmax), ylim=(-1.2*ymax, 1.2*ymax))
	plt.xlabel('Distance [AU]')
	plt.ylabel('Distance [AU]')
	startx = []
	starty = []
	colors = cm.rainbow(np.linspace(0, 1, numberOfBodies))
	scat = []
	ax.annotate('Time: %g yr' %(0.0), xy=(-xmax*1.3, -ymax*1.1))
	ax.annotate('%g yrs' %(time[-1]), xy=(xmax*1.25, -ymax*1.1))
	line, = ax.plot([], [], lw=3)
	for i in range(numberOfBodies):
		ax.annotate(bodies[i], xy=(xmax*0.9, ymax*(1-0.1*i)), color=colors[i])
		startx.append(i)
		starty.append(ymax*(1-0.1*i))
	scat = plt.scatter(startx, starty, s=200*dotSize, color=colors)

	# max frames when saving. otherwise buffer full hand program hangs
	maxFrames = 500
	frameDivider = timeSteps / maxFrames
	divider = timeSteps / numberOfFrames
	#print timeSteps, numberOfFrames

	# create different animators depending on save or not
	if (saveAnimation):
		anim = animation.FuncAnimation(fig, animatePointsSave, fargs=(frameDivider, fig, scat, line), frames=timeSteps/frameDivider, interval=1, blit=True, repeat=False)
	else:
		anim = animation.FuncAnimation(fig, animatePointsSave, fargs=(divider, fig, scat, line), frames=timeSteps/divider, interval=1, blit=True, repeat=False)
	
	# save animation or not
	if (saveAnimation):
		animationpath = groundpath + '/animations/'
		
		try:
			os.chdir(groundpath)
		except:
			print 'Failed to change directory to: %s' %(groundpath)
			sys.exit()
		try:
			os.chdir(animationpath)
		except:
			try:
				os.mkdir('animations')
				os.chdir(animationpath)
			except:
				'Failed to create subfolder /animations/. Please create it manually.'
				sys.exit()
		FFMpegWriter = animation.writers['mencoder']
		writer = FFMpegWriter(fps=30)
		animationFilename = filenameBase + '.mp4'
		print '\nSaving animation %s to \n%s' %(animationFilename, animationpath)
		anim.save(animationFilename, writer, extra_args=['-loglevel quiet'])
	else:
		plt.show()
	sys.exit()


# plotting position
plt.figure(0)
ax = plt.axes(xlim=(-1.6*xmax, 1.6*xmax), ylim=(-1.2*ymax, 1.2*ymax))
#ax = plt.axes(xlim=(-20, 20), ylim=(-10, 30))
start = 0
#if (x[0,0] == x[0, 1] and x[0,1] == x[0,2]):
#	start = 1
#	plt.plot([0], [0], 'o')
for i in range(start, numberOfBodies, 1):
	plt.plot(x[i,:], y[i,:])
#plt.axis('equal')
plt.xlabel('Distance [AU]')
plt.ylabel('Distance [AU]')
plt.legend(bodies)
savefilepos = filenameBase + '.pos.png'
plt.savefig(savefilepos)
plt.show()


if (False):
	plt.figure(0)
	ax = plt.axes(xlim=(-1.6*max(x[2,:]), 1.6*max(x[2,:])), ylim=(-1.2*max(y[2,:]), 1.2*max(y[2,:])))
	start = 0
	if (x[0,0] == x[0, 1] and x[0,1] == x[0,2]):
		start = 1
		plt.plot([0], [0], 'o')
	for i in range(start, numberOfBodies, 1):
		plt.plot(x[i,:], y[i,:])
	#plt.axis('equal')
	plt.xlabel('Distance [AU]')
	plt.ylabel('Distance [AU]')
	plt.legend(bodies)
	savefilepos = filenameBase + '.pos.png'
	plt.savefig(savefilepos)
	plt.show()


if (False):
	plt.figure(9)
	ax = plt.axes(xlim=(-1.6*max(x[5,:]), 1.6*max(x[5,:])), ylim=(-1.2*max(y[5,:]), 1.2*max(y[5,:])))
	start = 0
	#if (x[0,0] == x[0, 1] and x[0,1] == x[0,2]):
	#	start = 1
	#	plt.plot([0], [0], 'o')
	for i in range(start, 6, 1):
		plt.plot(x[i,:], y[i,:])
	#plt.axis('equal')
	plt.xlabel('Distance [AU]')
	plt.ylabel('Distance [AU]')
	plt.legend(bodies)
	savefilepos = filenameBase + '.inner.pos.png'
	plt.savefig(savefilepos)

	plt.figure(10)
	ax = plt.axes(xlim=(-1.6*max(x[0,:]), 1.6*max(x[0,:])), ylim=(-1.2*max(y[0,:]), 1.2*max(y[0,:])))
	plt.plot(x[0,:], y[0,:])
	plt.xlabel('Distance [AU]')
	plt.ylabel('Distance [AU]')
	plt.legend(['sun'])
	savefilepos = filenameBase + '.sun.pos.png'
	plt.savefig(savefilepos)


# plot velocities
plt.figure(1)
plt.plot(time, dx[1,:], time, dy[1,:], time, np.sqrt(dx[1,:]**2 + dy[1,:]**2))
plt.legend(['$v_x$', '$v_y$', '$|v|$'])
plt.xlabel('Time [yr]')
plt.ylabel('Velocity [AU/yr]')
savefilevel = filenameBase + '.vel.png'
#plt.savefig(savefilevel)

# plot forces
plt.figure(2)
plt.plot(time, Fx[1,:], time, Fy[1,:])
plt.legend(['$F_x$', '$F_y$'], loc='best')
plt.ylabel('Force')
plt.xlabel('Time [yr]')
savefileforce = filenameBase + '.force.png'
#plt.show()
#plt.savefig(savefileforce)


# plot energies
plt.figure(3)
plt.plot(time, Ek, time, Ep, time, Etot)
plt.legend(['$E_k$', '$E_p$', '$E_{tot}$'], loc='best')
plt.xlabel('Time [yr]')
plt.ylabel('Energy')
savefileenergy = filenameBase + '.energy.png'
#plt.savefig(savefileenergy)

# plot angular momentum
plt.figure(4)
plt.plot(time, Lx, time, Ly, time, Lz)
plt.legend(['$L_x$', '$L_y$', '$L_z$'], loc='best')
plt.xlabel('Time [yr]')
plt.ylabel('Angular momentum')
savefileang = filenameBase + '.ang.png'
#plt.savefig(savefileang)


#plt.show()
