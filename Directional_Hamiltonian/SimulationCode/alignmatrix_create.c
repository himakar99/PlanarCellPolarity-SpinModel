// This function creates the alignment matrix and stores the alignment data before mcmc simulation

/* Note:
First column is fraction of internally aligned cells
Second column is fraction of externally aligned cells
Third column is fraction of completely aligned cells
Fourth column denotes the direction of completely-aligned cells
Size: (dim*dim)X4
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

unsigned long long **alignmatrix_create(int **spin, int **nbr, int totalcell, FILE *fp3, FILE *fp4, FILE *fp5){
    
    int i;
    double c_int = 0.0;
    double c_ext = 0.0;
    double c_comp = 0.0;
    
    fprintf(fp3, "internal\texternal\tcomplete\n");
    
    unsigned long long **alignMatrix; // Matrix which stores data about internal, external, and complete alignment
    
    alignMatrix = malloc(sizeof(unsigned long long*)*totalcell);

    for (i=0; i<totalcell; i++){
        
        alignMatrix[i] =  malloc(sizeof(unsigned long long*)*4);
        alignMatrix[i][0] = 0;
        alignMatrix[i][1] = 0;
        alignMatrix[i][2] = 0;
        alignMatrix[i][3] = 9;
        
        //Checking internal alignment
        if(spin[i][0]*spin[i][1]+spin[i][1]*spin[i][2]+spin[i][2]*spin[i][3]+spin[i][3]*spin[i][4]+spin[i][4]*spin[i][5]+spin[i][5]*spin[i][0] == 2){
                alignMatrix[i][0] = 1;
            }
        
        //Checking external alignment
        if(spin[i][0]*spin[nbr[i][0]][3] + spin[i][1]*spin[nbr[i][1]][4] + spin[i][2]*spin[nbr[i][2]][5] + spin[i][3]*spin[nbr[i][3]][0] + spin[i][4]*spin[nbr[i][4]][1] + spin[i][5]*spin[nbr[i][5]][2] == -4){
            alignMatrix[i][1] = 1;
            }
        
        alignMatrix[i][2] = (alignMatrix[i][0] + alignMatrix[i][1])/2.0;
        
        if(alignMatrix[i][2] == 1){
            c_comp = c_comp + alignMatrix[i][2];
            
            //Checking directionality
            if(spin[i][0]+ 2*spin[i][1] == 1){
                alignMatrix[i][3] = 0; //Aligned upward
            }
            else if(spin[i][0]+ 2*spin[i][1] == 3){
                alignMatrix[i][3] = 1; //Aligned at 60 degrees from vertical
            }
            
            else if(spin[i][0]+ 2*spin[i][1] == 2){
                alignMatrix[i][3] = 2; //Aligned at 120 degrees from vertical
            }
            
            else if(spin[i][0]+ 2*spin[i][1] == -1){
                alignMatrix[i][3] = 3; //Aligned at 180 degrees from vertical
            }
            
            else if(spin[i][0]+ 2*spin[i][1] == -3){
                alignMatrix[i][3] = 4; //Aligned at 240 degrees from vertical
            }
            
            else if(spin[i][0]+ 2*spin[i][1] == -2){
                alignMatrix[i][3] = 5; //Aligned at 300 degrees from vertical
            }
        }
        else{
            alignMatrix[i][2] = 0;
            alignMatrix[i][3] = 9;
        }
        // TO BE COMMENTED FOR RUNNING ANALYSIS
        fprintf(fp4, "%llu\t", alignMatrix[i][2]);
        fprintf(fp5, "%llu\t", alignMatrix[i][3]);
         
    }
    // TO BE COMMENTED FOR RUNNING ANALYSIS
    fprintf(fp3, "%lf\t%lf\t%lf\n", c_int/totalcell, c_ext/totalcell, c_comp/totalcell);
    
    return alignMatrix;
}


