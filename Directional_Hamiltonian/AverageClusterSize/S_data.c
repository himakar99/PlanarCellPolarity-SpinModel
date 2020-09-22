//  This programme outputs the Average CLuster Size and it's standard error by combining data across multiple runs. For a particular value of m3 and DIM, it iterates through all m1-m2 combinations as well as combines data across multiple runs for the particular m1-m2 combination
//  The programme collects the counts of each cluster size (1 to L*L) iteratively, and then calculates the Ns and the subsequent standard errors iteratively as well.s
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
#define MAX(a,b) (((a)>(b))?(a):(b))

//--------------------------------INPUT STARTS--------------------------------//
#define DIM 512
#define RUNS 6 // No. of runs
#define DATA_SAMPLES 90
#define R_I 1 //First m1-m2 folder
#define R_F 1 //Final m1-m2 folder
char m3[100] = "2p2";
int run[RUNS] = {1,2,3,4,5,6}; // An array to allow selection of specific runs
//---------------------------------INPUT ENDS---------------------------------//

double calculateSD(long data[], double mean, int M);

int main(void){
    
    double time_spent = 0.0; //Variable describing time taken for execution of programme
    clock_t begin= clock();
    
    long totalcell = DIM*DIM;
    int M = DATA_SAMPLES*RUNS;
    
    char filepath[200]; // A dynamic variable used to update paths and open files
    char folder[200]; // Variable that will contain the number of m1-m2 folder in the form of a string
    char output_name[200]; // Variable that holds the name of the output file
    
    
    for(int f = R_I; f<= R_F; f++){
        
        
        FILE **inputs; // A double pointer of FILE type in order to open data of multiple runs and compile them
        inputs = calloc(sizeof(FILE*), RUNS);
        
        long *clusterCount;// Array that holds the sum of number of clusters for each size (1 to L*L) over M samples
        clusterCount = calloc(sizeof(long*), totalcell+1);
        
        double *clusterCountError;// Array that holds standard error of number of cluster of each size
        clusterCountError = calloc(sizeof(double*), totalcell+1);
        
        double *Ns; // Array that holds the Cluster Number Density distribution for cluster sizes 1 to L*L
        Ns = calloc(sizeof(double*), totalcell+1);
        
        double *Ns_error; // Array that holds the standard error of Cluster Number Density for cluster sizes 1 to L*L
        Ns_error = calloc(sizeof(double*), totalcell+1);
        
        //Opening and naming output file. Output format is S_(m3 value)_(m1-m2 folder no.).csv
        getcwd(output_name, sizeof(filepath));
        strcat(output_name, "/S_l");
        sprintf(folder, "%d", DIM);
        strcat(output_name, folder);
        strcat(output_name, "_");
        strcat(output_name, m3);
        strcat(output_name, "_");
        sprintf(folder, "%d", f);
        strcat(output_name, folder);
        strcat(output_name, ".csv");
        FILE *fpt;
        fpt = fopen(output_name, "w");
        
        long sum_NsS = 0;
        long sum_NsS2 = 0;
        
        double NsS_error = 0;
        double NsS2_error = 0;
        
        for(int r=0; r<RUNS; r++){
            // Opening the data files of multiple runs to combine them
            getcwd(filepath, sizeof(filepath));
            strcat(filepath, "/r");
            sprintf(folder, "%d", run[r]);
            strcat(filepath, folder);
            strcat(filepath, "/l");
            sprintf(folder, "%d", DIM);
            strcat(filepath, folder);
            strcat(filepath, "/" );
            strcat(filepath, m3);
            strcat(filepath, "/");
            sprintf(folder, "%d", f);
            strcat(filepath, folder);
            strcat(filepath, ".csv");
            inputs[r] = fopen(filepath, "r");
        }
        
        for(long s = 1; s<=totalcell; s++){
            int sample;
            
            long *clusterData; // An array that holds the counts of a particular cluster size over M samples.
            clusterData = calloc(sizeof(long*), M);
            
            for(int r=0; r<RUNS; r++){
                
                for(int a = 0; a<DATA_SAMPLES; a++){
                    sample = (r*DATA_SAMPLES) + a;
                    fscanf(inputs[r], "%li", &clusterData[sample]); //Reading the data
                    //printf("%li\t")
                    
                    if(clusterData[sample] > 0){
                        clusterCount[s] = clusterCount[s] + clusterData[sample]; //Summing to calculate Ns
                    }
                }
            }
            double mean;
            mean = (double)clusterCount[s]/M;
            clusterCountError[s] = calculateSD(clusterData, mean, M); // Calculating the standard deviation
            free(clusterData);
            clusterCountError[s] = clusterCountError[s]/sqrt(M); // Calculating the standard error
            
            Ns[s] = (double)clusterCount[s]/(totalcell*M);
            
            
            if(clusterCount[s] != 0){
                Ns_error[s] = Ns[s]*(clusterCountError[s]/clusterCount[s]);
            }else{
                Ns_error[s] = 0;
            }
            
            
            NsS_error = NsS_error + pow(Ns_error[s]*s, 2); // Numerator term for NsS
            NsS2_error = NsS2_error + pow(Ns_error[s]*s*s, 2); // Numerator term for NsS2
            
            sum_NsS = sum_NsS + clusterCount[s]*s;
            sum_NsS2 = sum_NsS2 + clusterCount[s]*s*s;
        }
        
        for(int r=0; r<RUNS; r++){
            fclose(inputs[r]); //Closing all the files opened across multiple runs for a particular m3 and m1-m2 combination
        }
        
        
        double Avg_cluster_size;
        double Avg_cluster_size_error;
        if(sum_NsS == 0){
            Avg_cluster_size = 0;
            Avg_cluster_size_error = 0;
        }else{
            Avg_cluster_size = sum_NsS2/sum_NsS;
            
            Avg_cluster_size_error = Avg_cluster_size*(sqrt(((NsS2_error)/pow(sum_NsS2,2)) + ((NsS_error)/pow(sum_NsS,2))));
        }
        
        fprintf(fpt, "%lf\t%lf", Avg_cluster_size, Avg_cluster_size_error);
        free(clusterCount);
        free(clusterCountError);
        free(Ns);
        free(Ns_error);
        fclose(fpt);
    
    }
    
    clock_t end = clock();
    time_spent += (double)(end-begin)/CLOCKS_PER_SEC;
    printf("Total time spent (sec): %f\n", time_spent);
    
    return 0;
}


// Utility function used to calculate SD of values of an array
double calculateSD(long data[], double mean, int M){
        double sum = 0.0, SD = 0.0;
        int i;

    for (i = 0; i<M; ++i){
        SD += pow(data[i] - mean, 2);
    }
    return sqrt(SD / M);
}
