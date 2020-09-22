This set of code is C implementation of directional hamiltonian problem based on Kamleswar's work.

It generates 6 output files:
1. a logfile with parameter information. The parameters are m1, m2, dim, randSpin

2. Energy_data: energy of the system PER cell at every sampling step. Also adds the initial energy per cell ate the very beginning.
Data is stored as ((mcmcSteps-burningPhase)/sampleInt)X1

3. Spin_data: spin of each edge of all the cells at the end of the simulation.
Data is stored as (dim*dim)X6

4. alignment_ratio: It stores fraction information of cell alignment per sampling step including the initial time point. Here definitions of alignments are as per Kamleswar’s paper. 
Data is stored as ((mcmcSteps-burningPhase)/sampleInt)X3

5. Complete_alignment: It stores which cells in the lattice are completely aligned and not at each of the sampling step. 0 = unaligned ; 1 = completely aligned. Here definition of complete alignment is as per Kamleswar’s paper.
Data is stored as ((mcmcSteps-burningPhase)/sampleInt)X1. Each row is a string containing the completeness of the alignment (0 or 1) of all cells, hence, length of each string is (dim*dim).

6. Aligment_direction: It stores alignment of each cell in the lattice at every sampling step. Directions are marked from 0 to 5 as decided by us. 9 is for unaligned cell.
Data is stored as ((mcmcSteps-burningPhase)/sampleInt)X1. Each row is a string containing the direction of alignment of all cells, hence, length of each string is (dim*dim).


