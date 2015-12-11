/***********************************************************************************/
/* Program: mortality.c
   By: Brad Duthie                                             
   Description: Kills off individuals from both the ID and OFF matrices
   Compile: gcc mortality.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void mortality(double **ID, double **OFF, int Liv, int l, int M, double Beta1){

    int i;
    double survival;

    /* First kill the older generation individuals in Liv */
    for(i=0; i<Liv; i++){ /* If individuals are too old (over M -- usually zero), die */
        if(ID[i][4]>=(M+1)){ /* They now have no location and -1 age */
            ID[i][2] = -1;
            ID[i][3] = -1;
            ID[i][4] = -1; 
        }
        if(ID[i][2] < 0){    /* And, for added measure, ensure they stay dead below */
            ID[i][4] = -1;
        }
    }    

    for(i=0; i<l; i++){ /* Next figure out if anyone dies in the new generation, OFF */
        survival   = exp(-1 * Beta1 * OFF[i][7]);
        OFF[i][2]  = survival; /* Put in column 2 for now (est. inbr. dep in sumstat)*/
        if(randunif()>survival){ /* Check to see if the individual dies */
            OFF[i][4] = -1;  /* If so, then add the -1 in column 4 */
        }
    } /* All juveniles should now have either survived, or died as juveniles */
}


