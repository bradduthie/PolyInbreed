/***********************************************************************************/
/* Program: babies.c
   By: Brad Duthie                                             
   Description: Uses the tables ID and Roff to make new offspring.
   Compile: gcc babies.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"

int babies(double **ID, int Clu, int i){

    int Off;

    Off  = 0;

    if(ID[i][8]<0){    /* If female did not find a mate  */
        Off = 0;       /* She produces no offspring      */
    }else{             /* If she did find a mate         */
        Off = Clu;     /* She produces Clu offspring     */
    }

    return Off;
}
