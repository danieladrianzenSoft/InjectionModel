

InjectionModel contains MATLAB code for a research-driven needle injection model consisting of continuum mechanics (CM) and mass transport (MT) theories

main.m is the main file that handles parameter initialization, execution of all physics, and plotting of results.

getDilatationLaplace.m, getDisplacementLaplace.m, getVelocityLaplace.m and getVelocityGradLaplace.m are functions that get the solution to the continuum mechanics theories described by Netti et al., 2003 in Laplace space. I later get the inverse Laplace of these functions in main.m to get the resulting dilatation, solid displacement, fluid velocity and fluid velocity gradient respectively. 

getConcetration_v2.m computes the concentration of injectate using the convection-diffusion equation and finite differences

Daniel Adrianzen