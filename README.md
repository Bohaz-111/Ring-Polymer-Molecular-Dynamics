# Ring-Polymer-Molecular-Dynamics (RPMD) — Fortran Toolkit

Comparison between **quantum mechanical**, **ring-polymer molecular dynamics**, and **classical** simulations of a particle in a simple 1D anharmonic potential. Implements RPMD following Craig & Manolopoulos (2004) with utilities to compute **Kubo-transformed real-time correlation functions**.

---

## Features

* Classical MD baseline (**n = 1** bead) and full **RPMD** (**n ≥ 1**)
* Kubo-transformed correlation functions (\tilde C_{AB}(t)) for position-dependent observables
* Correct short-time quantum limit and detailed-balance symmetry
* Modular integrators (velocity-Verlet / split-operator) and optional thermostats
* Reproducible examples for harmonic and anharmonic test systems

---

## Theory (brief)

### Quantum–classical isomorphism of the partition function

For a 1D system with (\hat H = \hat p^2/2m + V(\hat x)), the quantum canonical partition function
[
Z = \mathrm{tr}!\left[e^{-\beta \hat H}\right]
]
is isomorphic to the classical partition function of an (n)-bead ring polymer:
[
Z = \lim_{n\to\infty} Z_n,\qquad
Z_n = \frac{1}{(2\pi\hbar)^n}\
