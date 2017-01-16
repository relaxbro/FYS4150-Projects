FYS4150-project3 The Solar System
================

main.cpp
---
Set up solar system. Uses initial data from intitialdata/.
	
	Arguments
	None

Dependencies
- armadillo
- boost

Classes
- src/solarsystem.h
- src/solarsystemobject.h
- src/verlet.h
- src/rungekutta4.h

plotting.py
---
Plot or animate the orbits of all celestial objects. Can also plot velocities, forces, energies and angular momentum.

	Arguments
	1. filenameBase (base for data filenames, same as in main.cpp)
	2. animate (optional, use 'animate' if animation is wanted)
	3a. save (optional, use 'save' for saving animation, saves 500 frames 
	(maximumm becuase of buffer size))
	or
	3b. numberOfFrames (optional, int for number of frames displayed in animation)

convertToAU.py
---
Converts SI-units to AU, yr and solar masses. Initial values in initialdata/solarsystemSI.txt.
Writes multiple files to initialdata/ for each system.

	Arguments
	none
