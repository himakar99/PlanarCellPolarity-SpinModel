// This function generates the 2-dimensional matrix "spin", represented as "spin" in the program, of dimensions (dim^2 X 6)
// This is used to seed the intial spin cofigeration in the matrix. 
// seedSpin_func.c
//  seedSpin
//
//  Created by Himakar Sreerangam on 12/10/19.
//  Copyright © 2019 Himakar Sreerangam. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h> // GNU-Scientific Library introduced for random number generation.

// A function to generate a random permutation of arr[]. Uses Fisher-Yates Algorithm
void randomize1 (int arr[], int n, gsl_rng * r);

// A utility function to swap two integers
void swap1 (int *a, int *b);

int **seedSpin_func(int dim, int randSpin, unsigned long int mySeed, gsl_rng * r){
    
    int L = dim;
    int x;
    
    unsigned long int selectRow;
    
    selectRow = gsl_rng_uniform_int(r, 6);

    int **spin; // Integer double pointer that is passed back to main
    
    int allConfig[6][6] = {{-1, -1, 0, 1, 1, 0},
        {0, -1, -1, 0, 1, 1},
        {1, 0, -1, -1, 0, 1},
        {1, 1, 0, -1, -1, 0},
        {0, 1, 1, 0, -1, -1},
        {-1, 0, 1, 1, 0, -1}};
    
    spin = malloc(sizeof(int*)*(L*L)); // As dimensions of both matrix "nb" and "spin" are the same, we use the
    // same strategy as used in program "neb_func.c"
    
    for(x=0; x<(L*L); x++){
        spin[x] = malloc(sizeof(int*)*6);
        for (int i=0; i<6; i++){
            spin[x][i] = allConfig[selectRow][i]; // Each edge of the cell is given a spin.
        }
        if(randSpin == 1){             // If user wishes to have a completely random "seeded" state
            randomize1(spin[x], 6, r);     // Function called to shuffle a 1D array.
        }
    }
    return spin;
}

// Utility functions
void randomize1 (int arr[], int n, gsl_rng * r){
    
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = gsl_rng_uniform_int(r, i+1);
        
        // Swap arr[i] with the element at random index
        swap1(&arr[i], &arr[j]);
    }
}

void swap1 (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
