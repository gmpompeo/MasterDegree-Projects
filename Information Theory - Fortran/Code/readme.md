## Brief guide to code files

Here we present a list with a brief description of the purpose of each code file - further details can be found in the documentation within the code itself or in the .pdf report.

+ **initialization.f90**: file to initialize the adjacency matrix for each problem

+ **hamiltonians.f90**: file to build the hamiltonians (quantum and classical) of each NP problem

+ **energy.f90**: file to compute the "theoretical" ground state energy for the different algorithms

+ **quantum.f90**: file with implementation of Adiabatic Quantum Optimization algorithm (AQO)

+ **sa.f90**: file with implementation of classical simulated annealing (SA)

+ **sqa.f90**: file with implementation of simulated quantum annealing (SQA)
 
+ **checkpoint.f90**: file with checkpoint tools 

+ **debug.f90**: file with debugging subroutines (especially for AQO)

+ **utils.f90**: file with mathematical and general-purpose utilities 

+ **im.f90**: merging file with full implementation of Ising model (IM)

+ **gp.f90**: merging file with full implementation of graph partitioning problem (GP)

+ **vc.f90**: merging file with full implementation of vertex cover problem (VC)

+ **ts.f90**: merging file with full implementation of traveling salesman problem (TS)

+ **main.f90**: main program of the analysis

+ **analysis.py**: Python script working as launch pad for the analysis

***


The full analysis can be run by typing 

```python analysis.py``` 

on terminal. The script will ask the end-user a series of _y/n_ questions:
1) if he/she wants the Fortran analysis to be executed (the user can then pick what problem to run and this will prompt further problem-specific questions and parameter initializations);
2) if he/she wants the results to be shown for each problem, which will result in plots production.
Needed directories are created automatically.


