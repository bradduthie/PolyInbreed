/***********************************************************************************/
/* Program: sumstats.c
   By: Brad Duthie                                             
   Description: Returns summary stats for the inbreed program.
   Compile: gcc mortality.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"

void sumstats(double *RES, double **OFF, int M, int Nloci, int loadstart, int load, 
    int Neutstart, int Neutral, int Active, int l, int Imm, int xlen,
    double **REGR, int generation, int prevOff){

    int i, j;    
    double g, k, h, m, b, d;
    double sumxy, sumx, sumy, sumx2;    


    /* ==========================================================*/
    /* Calculate mean allele values     =========================*/
    /* ==========================================================*/        
    g    = 0;
    h    = 0;
    d    = 0;
    k    = 0;
    for(i=0; i<l; i++){
        if(OFF[i][4] >= 0 && OFF[i][4] <= M){
            for(j=0; j<Active; j++){
                g += OFF[i][((4*j)+10)];
                g += OFF[i][((4*j)+11)];
            }
            for(j=0; j<Neutral; j++){
                h += OFF[i][((4*j)+Neutstart)];
                h += OFF[i][((4*j)+Neutstart+1)];
            }
            for(j=0; j<load; j++){
                d += OFF[i][((4*j)+loadstart)];
                d += OFF[i][((4*j)+loadstart+1)];
            }
            k += 2;
        }
    }

    RES[0] = g / (Active * k);  /* Strategy allele frequency */
    RES[1] = h / (Neutral * k); /* Neutral allele frequency */
    RES[2] = d / (load * k);    /* Load allele frequency */


    /* ==========================================================*/
    /* Calculate Standard deviation values  =====================*/
    /* ==========================================================*/        
    g    = 0;
    h    = 0;
    d    = 0;
    k    = 0;
    for(i=0; i<l; i++){
        if(OFF[i][4] >= 0 && OFF[i][4] <= M){
            for(j=0; j<Active; j++){
                g  += (OFF[i][((4*j)+10)]-RES[0])*(OFF[i][((4*j)+10)]-RES[0]);
                g  += (OFF[i][((4*j)+11)]-RES[0])*(OFF[i][((4*j)+11)]-RES[0]);
            }
            for(j=0; j<Neutral; j++){
                h  += (OFF[i][((4*j)+Neutstart)]-RES[1])   *
                      (OFF[i][((4*j)+Neutstart)]-RES[1]);
                h  += (OFF[i][((4*j)+Neutstart+1)]-RES[1]) *
                      (OFF[i][((4*j)+Neutstart+1)]-RES[1]);
            }
            for(j=0; j<load; j++){
                d  += (OFF[i][((4*j)+loadstart)]-RES[2])   *
                      (OFF[i][((4*j)+loadstart)]-RES[2]);
                d  += (OFF[i][((4*j)+loadstart+1)]-RES[2]) *
                      (OFF[i][((4*j)+loadstart+1)]-RES[2]);
            }
            k += 2;
        }
    }

    RES[3] = sqrt((1/(Active  * k)) * g);  /* Strategy allele stdev */
    RES[4] = sqrt((1/(Neutral * k)) * h);  /* Neutral allele stdev */
    RES[5] = sqrt((1/(load    * k)) * d);  /* Load allele stdev */

    /* ==========================================================*/
    /* Regression coefficients -- juvenile survival =============*/
    /* ==========================================================*/    

    if(generation > 0){
        sumx   = 0; /* Sum of inbreeding coefficients */
        sumy   = 0; /* Sum of probability of survival */
        sumxy  = 0; /* Sum of inbreeding coefficient times prob survival */

        for(i=0; i<prevOff; i++){
            sumx   += REGR[i][1]; /* Individual's f coefficient */
            sumy   += log(REGR[i][6]); /* LOG probability of juvenile survival */
            sumxy  += REGR[i][1] * log(REGR[i][6]);
            sumx2  += REGR[i][1] * REGR[i][1];
        }

        /* Computing the slope below */
        m = (prevOff * sumxy - sumx * sumy) / (prevOff * sumx2 - sumx * sumx);
        /* Computing the intercept below */
        b = (sumy * sumx2 - sumx * sumxy) / (prevOff * sumx2 - sumx * sumx);
    
        RES[6]  = sumx / prevOff; /* Mean inbreeding coefficient */
        RES[7]  = sumy / prevOff; /* Mean log probability of juvenile survival */
        RES[8]  = m; /* Slope of juvenile inbreeding depression regression */
        RES[9]  = b; /* Intercept of juvenile inbreeding depression regression */
    }


}






