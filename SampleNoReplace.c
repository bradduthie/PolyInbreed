/***************************************************************/
/* Program: randunif.c
   By: Brad Duthie                                             
   Description: Will sample size "samp" from Vec without replace
   Compile: gcc randpois.c -ansi -Wall -pedantic    */
/***************************************************************/

#include "randunif.h" /* Need random unform numbers */
#include "array.h" /* Array used to handle vectors */

void SampleWithoutReplacement(int Pop, int samp, int *Vec){

    int n = samp;
    int N = Pop;

    int t = 0; /* total input records dealt with */
    int m = 0; /* number of items selected so far */    
    double u;

    while (m < n){
        u = randunif();
        if ((N - t)*u >= n - m ){
            t++;
        }else{
           Vec[m] = t;
           t++; 
           m++;
        }
    }
}

void shufflesamp(int *Vec, int samp){
    int i, j, temp;
    for (i = samp-1; i > 0; i--){
        j      = floor(randunif()*(i+1));
        temp   = Vec[i];
        Vec[i] = Vec[j];
        Vec[j] = temp;
    }
}


