//  This programme outputs the number of clusters of size s for each s over the number of data samples in each run. In order to calculate this it reads the file "alignment_direction.txt" for the run.
//  Created by Himakar Sreerangam on 02/10/19.
//  Copyright © 2019 Himakar Sreerangam. All rights reserved.
//
#include <stdbool.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nbr_create.c"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


//--------------------------------INPUT STARTS--------------------------------//
#define DIM 512
#define DATA_SAMPLES 90
#define R_I 2
#define R_F 3
//---------------------------------INPUT ENDS---------------------------------//

long countClusters(int M[], int *visited, long cluster_count[], long **clusterArray);

void DFS(int M[], long n_cell, int *visited, long count, long cluster_count[]);

int isSafe(int M[], long n_cell, int *visited);

void reverse(char *s);

void push(int item);

int pop();

int peek();

bool isStackEmpty();

// Creating nbr matrix
long **nbr;
int stack[DIM*DIM];
int top = -1;
int sample = 0;

int main(void){
    
    double time_spent = 0.0; //Variable describing time taken for execution of programme
    clock_t begin= clock();
    
    char filepath[200]; // A dynamic variable used to update paths and open files
    char folder[200]; // Variable that will contain the number of folder in the form of a strin
    
    for(int f = R_I; f<=R_F; f++ ){
        getcwd(filepath, sizeof(filepath));
        strcat(filepath, "/");
        sprintf(folder, "%d", f);
        strcat(filepath, folder);
        strcat(filepath, ".csv");
        
        //Opening output file
        FILE *fpt;
        fpt = fopen(filepath, "w");
        
        long totalcell = (DIM*DIM);
        nbr = nbr_create(DIM);
        
        //Printing column headers
       /* fprintf(fpt, "Sample\t");
        for(long a=1; a<=totalcell; a++){
            fprintf(fpt, "Size(%li)\t", a);
        }
        fprintf(fpt, "\n");*/
        
        getcwd(filepath, sizeof(filepath));
        strcat(filepath, "/");
        sprintf(folder, "%d", f);
        strcat(filepath, folder);
        strcat(filepath, "/");
        
        strcat(filepath, "alignment_direction.txt");
        FILE *fp1;
        printf("filename is %s\n",filepath);
        fp1 =  fopen(filepath, "r");
        int *num_initial; // Variable to extract alignment direction data at mcmcStep = 0;
        num_initial = calloc(sizeof(int*), totalcell);
            
        // Extracting initial data, i.e., mcmcStep = 0.
        for(long a = 0; a<DIM*DIM; a++){
            fscanf(fp1, "%d[^\t]", &num_initial[a]);
        }
        free(num_initial);
        // mcmcStep = 0 alignment direction data extracted
            
        char extract[totalcell+1]; //Variable used to extract sample alignment direction data
        long len;
        fscanf(fp1, "\t");
            
            
        sample = 0; // sample varilable that gets incremented
            
        int *num_array; // Variable used to convert sample alignment direction data from string to array
        num_array = calloc(sizeof(int*), totalcell);
        
        long **clusterArray;
        clusterArray = calloc(sizeof(long), totalcell+1);
            
        while(sample<DATA_SAMPLES && fgets(extract, totalcell+1, fp1)!= NULL){
                
            len = strlen(extract);
            if(len==totalcell){
                
                int *M;
                M = calloc(sizeof(int*), totalcell);
                int *visited;
                visited = calloc(sizeof(int*), totalcell);
                    
                for(long a = 0; a<len; a++){
                    num_array[a] = extract[a] - 48;
                    M[a] = 0;
                    visited[a] = 0;
                    if(sample == 0){
                        clusterArray[a] = calloc(sizeof(long), DATA_SAMPLES);
                    }
                    clusterArray[a][sample] = 0;
                    if(num_array[a] == 0){
                        M[a] = 1;
                    }
                }
                clusterArray[totalcell] = calloc(sizeof(long), DATA_SAMPLES);
                clusterArray[totalcell][sample] = 0;
                    
                long *cluster_count; //Array to keep track of the various clsuters in the lattice
                cluster_count = calloc(sizeof(long*), totalcell);
                    
                long max = 0;
                max = countClusters(M, visited, cluster_count, clusterArray); // Obtaining the size of every cluster in the lattice. Smallest cluster size = 1 (1 cell)
                
                for(long s = totalcell; s>0; s--){
                    if(s ==  max && clusterArray[s][sample] > 0){
                        clusterArray[s][sample] = 0;
                        break;
                    }
                }
                    
                long n_clusters = 0;
                int max_count = 0;
                free(M);
                free(visited);
                   
                sample++;
                free(cluster_count);

            } //If block
                
        } // While block
        fclose(fp1);
        free(num_array);
        
        for(long j = 1; j<=totalcell; j++){
            for (int i=0; i<DATA_SAMPLES; i++){
                fprintf(fpt, "%li\t", clusterArray[j][i]);
            }
            fprintf(fpt,"\n");
            free(clusterArray[j]);
        }
        free(clusterArray);
        fclose(fpt);
    } // For block
    
    
    clock_t end = clock();
    time_spent += (double)(end-begin)/CLOCKS_PER_SEC;
    printf("Total time spent (sec): %f\n", time_spent);
    
    return 0;
}

//Utility function to return an array containing clusters of different sizes
long countClusters(int M[], int *visited, long cluster_count[], long **clusterArray){
    
    // Initialize count as 0 and travese through the all cells of
    
    
    long count = 0;
    long max  = -1;
    for (long i = 0; i < DIM*DIM; i++){
        if (M[i] && !visited[i]){ // If a cell with value 1 is not // visited yet, then new cluster found
            cluster_count[count]++;
            DFS(M, i, visited, count, cluster_count); // Visit all cells in this cluster.
            max = MAX(max, cluster_count[count]);
            clusterArray[cluster_count[count]][sample]++;
            ++count; // and increment cluster count
        }
    }
    return max;
}

//Utility function to perform DFS for our nodes and edges
void DFS(int M[], long n_cell, int *visited, long count, long cluster_count[]){

    // Mark this cell as visited
    visited[n_cell] = 1;
    
    push(n_cell);

    while(!isStackEmpty()){
        
        n_cell = peek();
        pop();
        
        // Recur for all connected neighbours
        for (int k = 0; k < 6; ++k){
            if (isSafe(M, nbr[n_cell][k], visited)){
                cluster_count[count]++;
                visited[nbr[n_cell][k]] = 1;
                push(nbr[n_cell][k]);
            }
        }
    }
}

//Utility function to check if cell is 0-aligned, i.e., M[i] == 1 and not visited
int isSafe(int M[], long n_cell, int *visited){
// and not yet visited
    return (M[n_cell] && !visited[n_cell]);
}

//Utility function to reverse a string
void reverse(char *s){
   int length, c;
   char *begin, *end, temp;
 
   length = strlen(s);
   begin  = s;
   end    = s;
 
   for (c = 0; c < length - 1; c++)
      end++;
 
   for (c = 0; c < length/2; c++)
   {
      temp   = *end;
      *end   = *begin;
      *begin = temp;
 
      begin++;
      end--;
   }
}

//stack functions
void push(int item) {
   stack[++top] = item;
}

int pop() {
   return stack[top--];
}

int peek() {
   return stack[top];
}

bool isStackEmpty() {
   return top == -1;
}
