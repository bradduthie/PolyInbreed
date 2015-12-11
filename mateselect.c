/***********************************************************************************/
/* Program: mateselect.c
   By: Brad Duthie                                             
   Description: Allows the dominant sex to select a mate based on a set of criteria.
   Compile: gcc mateselect.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"

void mateselect(double **ID, double **Rmat, int Liv, int Nloci, int M, int mc,
    int Active, double alpha, int Kind, int wpe, double wpadj, int WpRestr){
    /* Note that ID should be sorted now, females at array top */    
    /* Note that females should be sorted by quality at top of ID */

    int i, j, k, h, ch, s, no, ccc, PrCount, g;
    double r, new, *PrM, *PrT, PrS, PrR, PrC;

    /* ===========================================================*/
    /* Determine where female and male sex row begins in ID array */
    /* ===========================================================*/        
    k = 0; /*Below counts the number of females */
    for(i=0; i<Liv; i++){
        if(ID[i][1]==0 && ID[i][4] > -1 && ID[i][4] <= M && ID[i][2] > -1){
            k++; 
        }else{
            break;
        }
    } /* So the females that can breed end at k */

    s = 0; /* Below gets the starting point for males */
    for(i=0; i<Liv; i++){
        if(ID[i][1]==0){
            s++;
        }else{
            break;
        }
    } /* So the males start at s */

    /* =======================================================================*/
    /* Each female individual assesses and chooses one or more mates          */
    /* =======================================================================*/

    MAKE_1ARRAY(PrM,Liv); /* Vector for the probability of selecting an individual */
    MAKE_1ARRAY(PrT,Liv); /* Cumulative probability vector used for sampling */

    for(i=0; i<k; i++){ /* For each female, starting with best */
        for(j=0; j<Liv; j++){ /* For each other individual in the population */
            PrM[j] = 0; /* Just start everything off with zero */
            PrT[j] = 0; /* The cumulative vector elements are zeros too */
        }
        ch  = -1; /* Initially, the individual has not made a mate choice(s) */
        no  = 0;  /* Used to break out of loop later if no more mates available */
        for(j=s; j<Liv; j++){ /* Check out each of opposite sex available */
            /* If potential mate not taken too many times or dead */
            if(ID[j][1]>(-1*mc) && ID[j][4]>-1 && 
            ID[j][4]<=M && ID[j][2] >= 0 && ID[j][3]>=0){
                if(Kind == 1){ /* If complete kin recognition */
                    r = Rmat[i][j+1]; /* Kinship with potential mate */
                }else{/* Below restricts to siblings */
                    if(ID[i][5]==ID[j][5] || ID[i][6]==ID[j][6]){
                        r = Rmat[i][j+1];
                    }else{ /* If not a sibling, treat as unrelated */
                        r = 0;
                    }
                }
                new = 0; /* Start with zero, but default will be one (below) */
                for(h=0; h<Active; h++){ /* Strategy alleles affect quality */
                    new += (r*alpha*ID[i][(h*4)+10]*wpe) + (r*wpadj); 
                    new += (r*alpha*ID[i][(h*4)+11]*wpe) + (r*wpadj); 
                } /* new is now higher if a relative is preferred, lower if disgust */
                if(new < 0){ /* If tends to inbreeding avoidance */
                    PrM[j] = 1 / (1 + (-1 * new));
                }else{ /* Now less likely to mate with rel than non-rel */
                    PrM[j] = 1 + new;
                } /* Else pos implies pref, and qual is just 1 + the val */
                no++; /* Adds recognition that a mate(s) was available */
            } /* Quality now increased (new) or decreased (1/new) with kinship */
        }
        if(no == 0){ /* Checks to see if a mate was available in that last loop */
            break; /* No need to continue if there are no mates left. */
        }
    
        /* ================================================================ */
        /* RESTRICT THE AVAILABILITY OF WP MATES TO A RANDOM SUBSET         */
        /* ================================================================ */

        if(WpRestr > 0){  /* Choice will be within restriction */
            ccc     = 100000; /* Just make sure there's no infinite loop */
            PrCount = 0; /* This will count possible mates (to be restricted) */
            for(j=0; j<Liv; j++){ /* Go through all of the individuals */
                if(PrM[j] > 0){ /* Check the probability this individual will sire */
                    PrCount++; /* Add one to PrCount if the probability is > 0 */
                } /* Now we have a count of all the possible sires */
            } /* Below, we will randomly restrict this count */
            while(PrCount > WpRestr){
                g = randunif() * Liv; /* Get a random sample */
                if(PrM[g] > 0){  /* If this was a possible sire */
                    PrM[g] = 0; /* Make him no longer a possible sire */
                    PrCount--;   /* And take one off of the counter */
                } /* At the end of the while loop, we should have only EpRestr */
                ccc--; /* Take one off of the possible 10k */
                if(ccc == 0){ /* If it got all the way through 10k iterations */
                    printf("\n\n INFINITE LOOP POSSIBLE \n\n");
                    break; /* Go through and check to make sure loop is not */
                } /* behaving ridiculously */
            } /* Males with the possibility of being chosen, comprising a restricted */
        } /* Subset that will be used for females to sample from */

        /* ================================================================ */
        /* SELECTION IS PROBABILISTIC BASED ON MATE QUALITY                 */
        /* ================================================================ */

        PrS = 0; /* Sum of the prob vector for choosing */
        for(j=0; j<Liv; j++){
            PrS += PrM[j];
        } /* PrS now can be thought of as `total quality' */
        if(PrS > 0){ /* If some probability of getting mate */
            PrR = 0; /* Running sum of probability vect */
            for(j=0; j<Liv; j++){
                PrT[j] = (PrM[j] / PrS) + PrR;
                PrR = PrT[j];
            } /* PrT now increasing vector, elements 0-1 */
            PrC = randunif(); /* Vector prob select*/
            j = 0; /* j increases to find position selected */
            while(PrT[j] < PrC){
                j++; /* Stops when random PrC pos found */
            }
            ch = j; /* The choice of individual is now j */
            PrS -= PrM[j];
            PrM[j] = 0;
            ID[i][8] = ch; /* Add i's choice into her col 8 */
        }else{ /* Else there is no probability of getting mate */
            ch = -1; /* So use -1 to indicate no mate acquired */
        }    
        if(ch > -1){ /* If i actually made some choice */
            ID[ch][8] = i; /* Change the mate */
            if(ID[ch][1] >= 0){ /* If mate not claimed  */
                ID[ch][1] = -1; /* Mate 1st time */
            }else{
                ID[ch][1]--; /* Else subtract one */
            } /* Above tracks times ch has been selected */
            ID[ch][2] = ID[i][2]; /* non to dom */    
            ID[ch][3] = ID[i][3]; /* non-dom to dom */
        }
    } 
    FREE_1ARRAY(PrM); /* Free the memory; vector is no longer needed */
    FREE_1ARRAY(PrT); /* Free the memory; vector is no longer needed */
}

































