FYS4150 project 2
================

main.cpp
---
Set up problem and use the Jacobi method to find eigenvalues.
	
	Arguments
	1. n (int size of matrix)
	2. numberOfElectrons (int number of electrons)
	3. rhoFinal (double rho max)
	4. omega (double omega frequenzy)

Uses
- src/eigenjacobirotation.cpp
- src/eigenjacobirotation.h

plotIterationsVSN.py
---
Plot number of iterations needed and calculation time.

	Arguments
	1. nMin (int starting value for size of matrix)
	2. nMax (int last value for n)
	3. rhoFinal (double rho max)

plotVSArma.py
---
Plot calculation time to compare with armadillo.

	Arguments
	1. nMin (int starting value for size of matrix)
	2. nMax (int last value for n)
	3. rhoFinal (double rho max)

plotOmega.py
---
Print and plot the first eigenvalue for multiple frequencies omega.

	Arguments
	1. n (int size of matrix)
	2. rhoFinal (double rho max)

plotWave.py
---
Plot the wavefunction for the three lowest eigenvalues. Uses the eigenvectors

	No arguments