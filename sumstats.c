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
    double g, k, h, d, x, y, Ex, Ey, Exy, sumx;

    /* ==========================================================*/
    /* Calculate mean allele values     =========================*/
    /* ==========================================================*/        
    g    = 0;
    h    = 0;
    d    = 0;
    k    = 0;
    for(i=0; i<l; i++){
        if(OFF[i][4] >= 0 && OFF[i][4] <= M){
            x = 0;
            y = 0;
            for(j=0; j<Active; j++){
                x += OFF[i][((4*j)+10)];
                x += OFF[i][((4*j)+11)];
            }
            g += x;
            for(j=0; j<Neutral; j++){
                y += OFF[i][((4*j)+Neutstart)];
                y += OFF[i][((4*j)+Neutstart+1)];
            }
            h += y;
            for(j=0; j<load; j++){
                d += OFF[i][((4*j)+loadstart)];
                d += OFF[i][((4*j)+loadstart+1)];
            }
            Ex  += x;
            Ey  += y;
            Exy += x*y;
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

        for(i=0; i<prevOff; i++){
            sumx   += REGR[i][1]; /* Individual's f coefficient */
        }
    
        RES[6]  = sumx / prevOff; /* Mean inbreeding coefficient */

    }

    /* Calculates the covariance between pre-cop inbreeding and polyandry */
    RES[7] = (Exy/(0.5*k)) - ( Ex/(0.5*k) * Ey/(0.5*k) ); /* Trait cov */
}






