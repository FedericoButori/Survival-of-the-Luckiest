# Survival-of-the-Luckiest
In this project I implemented in MATLAB a numerical scheme for the solution of a stochastic partial differential equation coming from population dynamics. Then I analyzed the effects of the random perturbations on the system

The system under study consists in an advection-diffusion-reaction system. The unknown is a function u(t,x) representing the density of some species at location x. 
The diffusion is modeled by D\Delta u(t,x), the advection by the usual transport term V(x)\cdot \nabla u(t,x) where V is some given vector field and the reaction by the logistic growth function r*u*(1-u/M) where M is the maximal capacity of the environment and r is the reproductive rate of the species. 
The reproductive rate is random and delta correlated in space, representing a random environment which can be fertile or hostile. 
The species is ultimately subject also to random advection caused by external perturbations (a strong wind or a turbulent environment). 
The goal is to understand how the random wind modifies the stationary equilibrium of the system, in particular analyzing the shape, elongation and size of the colonies which survive around fertile regions. 

# Files overlook

The repository consists in three matlab codes:
* A simple code solving a diffusion equation (the heat equation)
* A code solving advection-diffusion
* A final code solving a stochastic advection-diffusion-reaction with isotropic transport noise

Finally, a short pdf file describing in more details the motivation, analysis and results of the project, with additional figures.
