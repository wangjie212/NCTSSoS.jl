# [Polynomials in Non-Commutative Optimization](@id polynomials)

`NCTSSoS.jl` is designed to solve optimization problems where the variables do
not commute, a scenario that is the foundation of quantum mechanics. This page
outlines the key polynomial types used in the package and their specific
applications in many-body physics.

![Relation between different Polynomials](../assets/problem_classification.typ.svg)

## [Non-Commutative Polynomials: The Language of Quantum Systems](@id non-commutative-polynomial)

In quantum mechanics, the order of measurements matters. Non-commutative
polynomials provide the natural language for this reality. They are the building
blocks for constructing the **Hamiltonian** of a quantum system, which describes
its total energy.

## [Expectation Value Polynomials: Probing Non-Local Correlations](@id tracial-polynomial)

These polynomials represent the expected value of a combination of measurement
outcomes. A key application is the formulation of **linear Bell inequalities**
(e.g., the CHSH inequality). By finding the maximum value of these polynomials,
physicists can test the limits of quantum non-locality and entanglement, often
by assuming the system is in a maximally entangled state.

## [State Polynomials: Testing Advanced Quantum Phenomena](@id state-polynomial)

State polynomials are functions of the expectation values themselves. This
structure is essential for defining **nonlinear Bell inequalities**. These more
complex inequalities provide powerful and subtle tests of quantum mechanics,
allowing researchers to explore the boundary between the quantum and classical
worlds in greater detail.

## [Eigenvalue Polynomials: Finding the Energy Spectrum](@id eigenvalue-polynomial)

These polynomials are used to find the possible energy levels (eigenvalues) of a
quantum Hamiltonian. A primary application is to find the **ground state
energy**—the lowest possible energy of a many-body system. This is a fundamental
problem in condensed matter physics and quantum chemistry, with applications in
materials science and drug discovery.
