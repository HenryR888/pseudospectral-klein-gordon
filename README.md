# Pseudospectral Method for Breather Collisions in the Nonlinear Klein-Gordon Equation

This project implements a pseudospectral method to simulate breather solutions of the nonlinear Klein-Gordon equation:

$$
u_{tt} - u_{xx} + 4u - 2u^3 = 0
$$

The spatial derivative is computed using the Fast Fourier Transform, while the time evolution is performed using a centered finite difference scheme. The project first verifies the method on a single breather solution and then studies head-on collisions between two breathers for different half-separation distances.

## Project Summary

The main goals of this project are:

- to implement a pseudospectral method for the nonlinear Klein-Gordon equation;
- to simulate a single localised breather solution;
- to simulate head-on breather collisions;
- to compare the maximum collision amplitude for different half-separation distances.

The simulations show that the strongest interaction occurs when the breathers are initially superimposed. As the initial half-separation distance increases, the collision amplitude decreases and then stabilises.

## Numerical Method

The spatial domain is discretised on $$[-L,L]$$, with $$L = 100$$ and $$N = 512$$ grid points. The second spatial derivative is computed in Fourier space using:

$$
\widehat{u_{xx}} = -q^2 \hat{u}.
$$

The time evolution uses a centered finite difference scheme:

$$
u_j^{k+1} = 2u_j^k - u_j^{k-1} + \tau^2 \left( (u_{xx})_j^k - 4u_j^k + 2(u_j^k)^3 \right).
$$

Full details and results are documented in [Project 2](project/Project%202.pdf).

## Attribution

The problem was provided as coursework in:

**MAM3042F – Advanced Numerical Methods**  
**University of Cape Town**

Lecturer: **Dr N. V. Alexeeva**

All implementations, experiments, and analysis in this repository are my own.