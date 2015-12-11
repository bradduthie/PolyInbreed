/***********************************************************************************/
/* Program: inbrdep.c 
   By: Brad Duthie  
   
   Description: Two functions to transfer offspring information and check survival 
   Compile: gcc inbrdep.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

/* This function transfers information from OFF to REGR, as needed */
void RegrTr(double **REGR, double **OFF, int prevOff){

    int i;

    for(i=0; i<prevOff; i++){
        REGR[i][0] = OFF[i][0]; /* Get individual number */
        REGR[i][1] = OFF[i][7]; /* This will store an individual f coefficient */
        REGR[i][2] = 0; /* Empty now, use for something later? */
        REGR[i][3] = OFF[i][1];
        REGR[i][4] = 1; /* Start assuming complete survival probability, adj later */
        REGR[i][5] = 0; /* Assume dead unless proved otherwise by becoming an adult */
        REGR[i][6] = OFF[i][2]; /* This was the juvenile probability of survival */
        REGR[i][7] = 0; /* An individual's offspring production will go here */
    }
}

/* This function checks to see if individuals from original OFF ultimately survived */
void inbrdep(double **ID, double **REGR, int prevOff, int Liv, int M){

    int i, j;

    for(i=0; i<Liv; i++){
        if(ID[i][4] > -1 && ID[i][4] <= M){ /* If i is alive */
            for(j=0; j<prevOff; j++){ /* Go through REGR and find it, j */
                if(ID[i][0] == REGR[j][0]){
                    REGR[j][5] = 1; /* Add a one to indicate survival */
                    break; /* No need to keep through this loop */
                }
            }
        }
    }
}


/* Prints kinship relative to an indivual to a kinship file (use with pedigree) */
void printMeanKinship(double **Rmof, double **ID, int Liv,
    int Active, int Neutral, int Neutstart, int load, int loadstart, int pid,
    double Beta1, double Pcost, int WpRestr, int EpRestr){

    int    j, k, i;
    double a, b, c;
    FILE *kinshipmeans; /* This will print out if Pedigree is also printing */

    kinshipmeans = fopen("kinshipmeans.txt","a+");

    for(j=0; j<Liv; j++){
        if(ID[j][4]==0){
            a = 0;
            b = 0;
            c = 0;
            for(i=0; i<Active; i++){
                a += ID[j][((4*i)+10)];
                a += ID[j][((4*i)+11)];
            }
            for(i=0; i<Neutral; i++){
                b += ID[j][((4*i)+Neutstart)];
                b += ID[j][((4*i)+Neutstart+1)];
            }
            for(i=0; i<load; i++){
                c += ID[j][((4*i)+loadstart)];
                c += ID[j][((4*i)+loadstart+1)];
            }
            for(k=0; k<Liv; k++){
                if(ID[k][4]==0 && j != k){
                    fprintf(kinshipmeans,"%d\t%d\t",WpRestr,EpRestr);
                    fprintf(kinshipmeans,"%f\t%f\t",Pcost,Beta1);
                    fprintf(kinshipmeans,"%d\t",pid);
                    fprintf(kinshipmeans,"%f\t",ID[j][0]); /* ID */
                    fprintf(kinshipmeans,"%f\t",ID[j][1]); /* Sex */
                    fprintf(kinshipmeans,"%f\t",a); /* WP-Mating */
                    fprintf(kinshipmeans,"%f\t",b); /* Polyandry */
                    fprintf(kinshipmeans,"%f\t",c); /* EP-mating */
                    fprintf(kinshipmeans,"%f\t",ID[k][0]); /* ID */
                    fprintf(kinshipmeans,"%f\n",Rmof[j][k+1]); /* Kinship */   
                }
            }
        }
    }
    fclose(kinshipmeans);
}




