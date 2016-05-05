/***************************************************************/
/* Program: randpois.c
   By: Brad Duthie                                             
   Description: Will produce a Poisson distributed random number
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***************************************************************/


#include "randunif.h" /* Need randunif header file */

int randpois(double lambda){
    double L, u, p;
    int k;

    p = 1;
    k = 0;
    L = exp(-1*lambda);
    do{
        k++;
        u = randunif();
        p = p * u;
    }while(p > L);
    return k-1;
}
