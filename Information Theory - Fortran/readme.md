## Information Theory

The project developed in this course required applications of computational methods to physical concepts and models.

In particular, NP-problems (such as graph partitioning, traveling salesman and more) were described with an Ising formalism. The idea was to optimize these systems as if they were physical, i.e. configurations of spins whose ground state needs to be found.

For this purpose, three algorithms were implemented: *adiabatic quantum optimization*, which simulates a quantum computer in the attempt to reach ground state through an adiabatic transformation, *simulated annealing* in a Metropolis framework and *quantum simulated annealing*, employing Trotter replicas to simulate quantum tunneling effects. This last algorithm, combining the flexibility of MC methods and the power of quantum tricks, proved to be the most efficient across several scenarios.