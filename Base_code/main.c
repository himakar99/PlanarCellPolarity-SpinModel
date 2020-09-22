// 
//  main.c
//  Main function to perfomr MCMC simulations for PCP model
//  as perfomed in Kamlswar's paper.
//  Created by Himakar Sreerangam on 20/10/19.
//  Copyright © 2019 Himakar Sreerangam. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h> // GNU-Scientific Library introduced for random number generation.
#include <time.h>
#include <sys/time.h> // Library introduced to ensue unique seed every simulation
#include "seedSpin_func.c"
#include "neb_func.c"
#include "mcmc.c"
#include "alignmatrix_create.c"
#include "alignmatrix_func.c"
#include <sys/stat.h>
#include <sys/types.h>


// A utility function to generate a unique seed for random generators every simulation
unsigned long int random_seed();

// A utility function to calculate the energy of the lattice
double energy_calc(int **nbr, int **spin, int totalcell, double m1, double m2);

int main(){
    
    double time_spent = 0.0; //Variable describing time taken for execution of programme
    clock_t begin= clock();
    
    unsigned long int mySeed; // Seed to ensure randomness
    mySeed = random_seed();
    
    char filename[30000]; // Character array to name files
    struct tm *timenow;
    
    //---------------------------------------------------------------------INPUTS START---------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    // Interaction parameters. m1 positive & m2 negative
    double m1,m2; 
    m1 = 4;
    m2 = -4;

    // dimension of lattice
    unsigned long long dim = 16; 

    // Intitial spin configuration. 
    //0 = all cells identically completly algined; 1 = random
    int randSpin = 1; 

    //Total number of mcmc-steps; change numerical value ONLY
    unsigned long long mcmcSteps = (dim*dim)*(long long)100000; 

    // Number of initial MC steps to be discarded; change numerical value ONLY
    unsigned long long burningPhase = (dim*dim)*(long long)2000; 

    // Data sampling interval AFTER burning phase;change numerical value ONLY
    unsigned long long sampleInt = (dim*dim)*(long long)1000; 
    

    //DO NOT Change the following unless required. Decides number temporary files for data
    unsigned long long filecount = 10000000; 
    
    
    //---------------------------------------------------------------------Output files---------------------------------------------------------------------//
    FILE *fp1; //Opening a "txt" file to write energy-data
    fp1 = fopen("Energy_data.txt", "w");
    
    FILE *fp2; //Opening a "txt" file to write spin-data
    fp2 = fopen("Spin_data.txt", "w");
    
    FILE *fp3; //Opening a "txt" file to write alignment-ratio data
    fp3 = fopen("alignment_ratio.txt", "w");
    
    FILE *fp4; //Opening a "txt" file to write alignment data i.e. 0 = unaligned ; 1 = aligned
    fp4 = fopen("complete_alignment.txt", "w");
    
    FILE *fp5; //Opening a "txt" file to write alignment-direction data
    fp5 = fopen("alignment_direction.txt", "w");
    
    FILE *fp6; //Opening a "txt" file to write log file parameter information 
    fp6 = fopen("logfile.txt", "w");
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------INPUTS END---------------------------------------------------------------------//
    
    int totalcell = dim*dim; // Total number of cells in the lattice
    int **spin; // Double pointer used to access 2-dimensional array
    int **nbr; // Double pointer used to access 2-dimensional array
    unsigned long long **alignMatrix; // Triple pointer used to access 3-dimensional array
    
    fprintf(fp6, "dim = %llu\n", dim);
    fprintf(fp6, "randspin = %d\n", randSpin);
    fprintf(fp6, "m1 = %f\n", (m1));
    fprintf(fp6, "m2 = %f\n", (m2));
    fprintf(fp6, "Total Number of MC steps = %llu\n", mcmcSteps);
    fprintf(fp6, "Total Number of burning steps = %llu\n", burningPhase);
    fprintf(fp6, "Sampling interval = %llu\n", sampleInt);
    fprintf(fp6, "data samples = %llu\n", (mcmcSteps-burningPhase)/sampleInt);
    
    
    // Random number generation set-up
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_ranlxs0;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    // Random number generation set-up completed
    
    //Creating spin, nbr, alignment matrices
    spin = seedSpin_func(dim, randSpin, mySeed, r); // Creating spin matrix
    nbr = neb_func(dim); // Creating nbr matrix
    alignMatrix = alignmatrix_create(spin, nbr, totalcell, fp3, fp4, fp5); // Creating alignment matrix
    //Matrices created
    
    // Calculating total initial Energy
    double energy;
    energy = energy_calc(nbr, spin, totalcell, m1, m2);
    fprintf(fp1, "%lf\n", (energy/totalcell));
    //Initial energy calculated
    
    //MCMC
    unsigned long long i = 0;
    
    for(i = 1; i<=burningPhase; i++){
        energy = energy + mcmc(spin, m1, m2, dim, totalcell, nbr, r);
    }
    
    for(i=burningPhase + (long long)1; i<=mcmcSteps; i = i + (long long)1){
        
        energy = energy + mcmc(spin, m1, m2, dim, totalcell, nbr, r);
        
        if(i%sampleInt == 0){
            fprintf(fp1, "%lf\n", (energy/totalcell));//Energy data is collected every sampleInt time-steps
            alignmatrix_func(alignMatrix, spin, nbr, totalcell, fp3, fp4, fp5);
        }
    }
    //MCMC completed
    printf("MCMC Simulation is over\n");
    
    gsl_rng_free (r); //Random number generator closed
    fclose(fp1);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    
    // For loop to write the spin matrix-data into the "txt" file
    for(int i=0; i<dim*dim; i++){
        for(int j=0; j<6;j++){
            fprintf(fp2, "%d", spin[i][j]);
        }
        fprintf(fp2, "\n");
    }
    fclose(fp2); //spin-matrix data printed
    
    clock_t end = clock();
    
    time_spent += (double)(end-begin)/CLOCKS_PER_SEC;
    fprintf(fp6, "Total time spent (sec): %f\n", time_spent);
    fclose(fp6);
    return 0;
}
//Main of program ends

// Utility functions
unsigned long int random_seed(){
    // To generate a unique seed for random number generator
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
    
}

double energy_calc(int **nbr, int **spin, int totalcell, double m1, double m2){
    double energy = 0.0;
    double h2 = 0;
    double h1 = 0;
    
    for (int i=0; i<totalcell; i++){
        h2 = h2 + m2*(spin[i][0]*spin[i][1]+ spin[i][1]*spin[i][2]+ spin[i][2]*spin[i][3]+ spin[i][3]*spin[i][4]+ spin[i][4]*spin[i][5]+ spin[i][5]*spin[i][0]);
        
        for (int j =0; j<6; j++){
            h1 = h1 + m1*(spin[nbr[i][j]][(j+3)%6]*spin[i][j]);
        }
    }
    energy = h1/2.0 + h2;
    return energy;
}
