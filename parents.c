/***********************************************************************************/
/* Program: parents.c
   By: Brad Duthie                                             
   Description: Determines which offspring came from which parents, and properties.
   Compile: gcc babies.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "array.h"
#include "randunif.h"
#include "randpois.h"
#include "randnorm.h"

void parents(double **ID, double **OFF, double **Rmof, int *O, int Nloci, int Ind, 
    int k, int l, int Clu, int Liv, int gener, double Beta1, double alpha,
    double *RES, int rep, int loadstart, int load, int Active, int prP, int Neutstart,
    int Neutral, int Kind, int M, int poe, int epe, double poadj, double epadj, 
    double mu, double mumu, double musd, int conSt, int EpAvail, int pid, int msel,
    int lastgen, int EpRestr, double Pcost, int condk){

    int i, j, g, h, m, s, MomP, DadP, extra, nmale, ch, PrCount, ccc, socmat;
    int *MALES;
    double pscore, nscore, mscore, gPy, r, PrS, PrR, PrC, new, a, ksoc;
    double *PrM, *PrT, *Mse;
    FILE *Pedigree, *Mates;

    /*===========================================================================*/
    /* Social male in mating pool, need to select others?  ======================*/
    /*===========================================================================*/
    h     = 0;
    extra = 0; /* Check to see if we need to find a new extra pair mate for nest */
    MAKE_1ARRAY(MALES,2); /* Just so to not free something that doesn't exist */
    MAKE_1ARRAY(Mse,2); /* Same as above, will only use if EpAvail > 0 */
    for(i=0; i<l; i++){
        if(i==O[h]){   /* If offspring # hits count of off nest*/    
            h++;       /* Move to the next nest */
            extra = 0; /* Need to find a new extra pair mate */
        }
        /* -- If new extra-pair males need to be selected -- who are they? */
        if(extra == 0){ /* If a new males not yet selected as an epm */
            gPy = 0.0;  /* Start out with a zero value for polyandry-like */ 
            if(poe == 1){ /* no change (gPy = 0) if poe=0 */
                for(j=0; j<Neutral; j++){
                    gPy += ID[h][((4*j)+Neutstart)]; 
                    gPy += ID[h][((4*j)+Neutstart+1)];
                }
            }
            gPy += poadj; /* Adjusts polyandry externally (in PolyIn.c) */
            if(condk == 1){ /* If prob of polyandry is conditional upon k */
                socmat = ID[h][8]+1;
                ksoc   = Rmof[h][socmat]; /* Get the k of the initial mate mate */
                new = 0; /* This and below just checks if avoiding (<0) or prefering (>0) */
                for(g=0; g<Active; g++){ /* Strategy alleles affect quality */
                    new += ID[h][(4*g)+10]*epe + epadj; 
                    new += ID[h][(4*g)+11]*epe + epadj; 
                } /* Now if preferring & cond == 1, adjust condition for polyandry */
                if(new > 0){ /* If preferring, want females to engange in polyandry */
                    ksoc = 1 - ksoc; /* when mate is more distant relative */
                } /* So subtract here so that high kinship between mates reduces polyandry */
            }else{
                ksoc = 1.0;
            }
            if(gPy < 0){ /* Lambda value has to be zero or positive */
                gPy = 0.0;
            }
            nmale = 1 + randpois(gPy*ksoc); /* So we need an array with each male */
            if(nmale > k){ /* If females more polyandrous than there are males */
                nmale = k; /* Cap the polyandry at the number of paired males */
            } /* Now males cannot be selected additionally as extra-pair mates */
            if(EpAvail > 0 && (EpAvail+1) < nmale){ /* If random subset */
                nmale = EpAvail; /* But would be more polyandrous than subset size */
            } /* Must truncate polyandry so female can't take more males than subset */
            if(EpAvail <= 0 && EpRestr > 0 && (EpRestr+1) < nmale){ /* If subset restricted */
                nmale = EpRestr; /* But would be more polyandrous than the restriction */
            } /* Must truncate polyandry so females can't take more than restriction */
            if(gPy <= 0.0 && nmale != 1){
                printf("ERROR IN PARENTS.C -- too many males selected\n");
            }
            if(nmale > 1){ /* If there is more than one male to consider */
                MAKE_1ARRAY(PrM,Liv); /* Pr of selecting a male */
                MAKE_1ARRAY(PrT,Liv); /* Cumulative prob for sampling */
                FREE_1ARRAY(MALES);
                MAKE_1ARRAY(MALES,nmale); /* Here is where the males are stored */
                for(j=0; j<Liv; j++){ /* For each other individual in the population */
                    PrM[j] = 0; /* Just start everything off with zero */
                    PrT[j] = 0; /* The cumulative vector elements are zeros too */
                }
                s = 0; /* Need to identify and select other males */
                for(j=0; j<Liv; j++){
                    if(ID[j][1]==0){
                        s++;
                    }else{
                        break;
                    }
                } /* Have identified now that the males start at s */
                /* ----- Everything below checks out new mates -----------------  */
                for(j=s; j<Liv; j++){ /* Check out each of opposite sex available */
                    if(ID[j][4]>-1 && ID[j][4]<=M && ID[j][2] >= 0 && ID[j][3] >= 0){
                        if(EpAvail > 0){ /* If EpAvail is positive, using a subset */
                            PrM[j] = 1; /* And the subset should be random */
                        }else{ /* Else we choose nmales based on genotype */
                            if(Kind == 1){ /* If complete kin recognition */
                                r = Rmof[h][j+1]; /* Kinship with potential mate */
                            }else{/* Below restricts to siblings */
                                if(ID[h][5]==ID[j][5] || ID[h][6]==ID[j][6]){
                                    r = Rmof[h][j+1];
                                }else{ /* If not a sibling, treat as unrelated */
                                    r = 0;
                                }
                            }
                            new = 0; /* Start with zero, but default will be one (below) */
                            if(conSt == 0){
                                for(g=0; g<load; g++){ /* Strategy alleles for qual */
                                    new += (r*alpha*ID[h][((4*g)+loadstart)]*epe)+(r*epadj); 
                                    new += (r*alpha*ID[h][((4*g)+loadstart+1)]*epe)+(r*epadj); 
                                }
                            }else{ /* epadj adjusts ep preference externally (PolyIn.c) */
                                for(g=0; g<Active; g++){ /* Strategy alleles affect quality */
                                    new += (r*alpha*ID[h][(4*g)+10]*epe) + (r*epadj); 
                                    new += (r*alpha*ID[h][(4*g)+11]*epe) + (r*epadj); 
                                } /* Note: Still adjusting with epe and epadj, not wp adjs */
                            }
                            if(new < 0){ /* If tends to inbreeding avoidance */
                                PrM[j] = 1 / (1 + (-1 * new));
                            }else{ /* Now less likely to mate with rel than non-rel */
                                PrM[j] = 1 + new;
                            } /* Else pos implies pref, and qual is just 1 + the val */
                            /* If nmale choice restricted to some number, but */
                        } /* Quality now increased (new) or decreased (1/new) with kinship */
                    } /* Like in mateselect, now move to actual selection from prob */
                } /* Below: If not using a subset, but want to restrict access to males */
                if(EpAvail <= 0 && EpRestr > 0){  /* Choice will be within restriction */
                    ccc     = 100000; /* Just make sure there's no infinite loop */
                    PrCount = 0; /* This will count possible mates (to be restricted) */
                    for(j=0; j<Liv; j++){ /* Go through all of the individuals */
                        if(PrM[j] > 0){ /* Check the probability this individual will sire */
                            PrCount++; /* Add one to PrCount if the probability is > 0 */
                        } /* Now we have a count of all the possible sires */
                    } /* Below, we will randomly restrict this count */
                    while(PrCount > EpRestr){
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
                MALES[0] = ID[h][8]; /* Social mate gets the first position */
                g        = ID[h][8]; /* g assigned as integer for line below */ 
                PrM[g]   = 0;        /* Social mate cannot take another position */
                for(j=1; j<nmale; j++){ /* Every other position (below) is EPMs */
                    PrS = 0; /* Sum of the prob vector for choosing */
                    for(g=0; g<Liv; g++){
                        PrS += PrM[g];
                    } /* PrS now can be thought of as `total quality' */
                    if(PrS > 0){ /* If some probability of getting mate */
                        PrR = 0; /* Running sum of probability vect */
                        for(g=0; g<Liv; g++){
                            PrT[g] = (PrM[g] / PrS) + PrR;
                            PrR = PrT[g];
                        } /* PrT now increasing vector, elements 0-1 */
                        PrC = randunif(); /* Vector prob select*/
                        g = 0; /* j increases to find position selected */
                        while(PrT[g] < PrC){
                            g++; /* Stops when random PrC pos found */
                        }
                        MALES[j] = g; /* The choice of individual is now g */
                        PrM[g]   = 0; /* g cannot be selected again */
                    }else{
                        MALES[j] = ID[h][8]; /* Not enough EPM? */
                    } /* Should never happen, but default to social pair */
                } /* Note, males cannot be exhausted. If alive, they're fair game */
                /* ----- Now we have a vector of wpm and epms -----------------  */
                FREE_1ARRAY(PrM); /* Note: If EpAvail > 0, MALES is a random subset */
                FREE_1ARRAY(PrT); /* Selection within non-random if EpAvail>0 (see below */
                /* If EpAvail > 0, EP males are random, selection occurs within subset */
                if(EpAvail > 0){ /* If true, males will be social + random set */
                    MAKE_1ARRAY(PrM,nmale); /* Use now for within the subset */
                    FREE_1ARRAY(Mse); /* Should always be one to free */
                    MAKE_1ARRAY(Mse,nmale); /* Cumulative prob for sampling */
                    for(j=0; j<nmale; j++){
                        g = MALES[j]; /* MALE[j] is the row of the male, what we need */
                        if(Kind == 1){ /* If complete kin recognition */
                            r = Rmof[h][g+1]; /* Kinship with potential mate */
                        }else{/* Below restricts to siblings */
                            if(ID[h][5]==ID[g][5] || ID[h][6]==ID[g][6]){
                                r = Rmof[h][g+1];
                            }else{ /* If not a sibling, treat as unrelated */
                                r = 0;
                            }
                        } /* Note that g is being recycled below, no longer position */
                        new = 0; /* Start with zero, but default will be one (below) */
                        if(conSt == 0){
                            for(g=0; g<load; g++){ /* Strategy alleles for qual */
                                new += (r*alpha*ID[h][((4*g)+loadstart)]*epe)+(r*epadj); 
                                new += (r*alpha*ID[h][((4*g)+loadstart+1)]*epe)+(r*epadj); 
                            }
                        }else{ /* epadj adjusts ep preference externally (PolyIn.c) */
                            for(g=0; g<Active; g++){ /* Strategy alleles affect quality */
                                new += (r*alpha*ID[h][(4*g)+10]*epe) + (r*epadj); 
                                new += (r*alpha*ID[h][(4*g)+11]*epe) + (r*epadj); 
                            } /* Note: Still adjusting with epe and epadj, not wp adjs */
                        }
                        if(new < 0){ /* If tends to inbreeding avoidance */
                            PrM[j] = 1 / (1 + (-1 * new));
                        }else{ /* Now less likely to mate with rel than non-rel */
                            PrM[j] = 1 + new;
                        } /* Else pos implies pref, and qual is just 1 + the val */
                    } 
                    PrS = 0; /* Sum of the prob vector for choosing */
                    for(g=0; g<nmale; g++){
                        PrS += PrM[g];
                    } /* PrS now can be thought of as `total quality' */
                    if(PrS > 0){ /* If some probability of getting mate */
                        PrR = 0; /* Running sum of probability vect */
                        for(g=0; g<nmale; g++){
                            Mse[g] = (PrM[g] / PrS) + PrR;
                            PrR = Mse[g];
                        } /* PrT now increasing vector, elements 0-1 */
                    } /* Will use PrT vector later in selection */
                    FREE_1ARRAY(PrM);
                }
            } /* Now we also know which males can sire the offsring of i */
            extra++; /* Don't need a new set of males now */
            if(msel == 1 && gener == (lastgen - 1)){ /* If printing females mate choice */
                Mates = fopen("mates.txt","a+"); /* Open the mates file */
                fprintf(Mates,"%d\t%f\t%f\t%f\t",pid,Beta1,Pcost,ID[h][0]);/*Print pid,mumID */
                fprintf(Mates,"%f\t%d\t",gPy,nmale); /* Print off tendancy, nmales */
                if(nmale == 1){ /* If she only selected one male */
                    ch   = ID[h][8];
                    fprintf(Mates,"%f\t",ID[ch][0]); /* Print the mate ID */
                    DadP = 0; /* Need to sweep through to find Rmof position of male */
                    while(ID[ch][0]!=Rmof[DadP][0]){
                        DadP++; /* Keep increasing until find male position in Rmof */
                    } /* Now can print the kinship next to the ID of the potential mate */
                    fprintf(Mates,"%f\t",Rmof[h][DadP+1]); /* Print the kinship */
                }else{ /* Else if she took more than one mate */
                    for(j=0; j<nmale; j++){ /* For each male the mother selected as a mate */
                        fprintf(Mates,"%f\t",ID[MALES[j]][0]); /* Print the mate ID */
                        DadP = 0; /* Need to sweep through to find Rmof position of male */
                        while(ID[MALES[j]][0]!=Rmof[DadP][0]){
                            DadP++; /* Keep increasing until find male position in Rmof */
                        } /* Now can print the kinship next to the ID of the potential mate */
                        fprintf(Mates,"%f\t",Rmof[h][DadP+1]); /* Print the kinship */
                    } /* Note: First male position will be the initially selected male */
                } /* Lastly, skip to the next line so a new female will be printed */
                fprintf(Mates,"\n");
                fclose(Mates);
            } /* Should now have a ragged array in the file `Mates.txt' */
        } /* Move along to actually dispensing paternity */
        if(nmale==1){ /* If there is no extra-pair paternity, m just social male */
            m  = ID[h][8];
        }else{ /* Else, check if we are operating with EpAvail */
            if(EpAvail > 0){ /* If so, choose within subset based on quality */
               PrC = randunif(); /* Vector prob select*/
               ch  = 0; /* j increases to find position selected */
               while(Mse[ch] < PrC){
                   ch++; /* Stops when random PrC pos found */
               }
               m = MALES[ch]; /* The choice of individual is now m */
            }else{ /* Else take a random subset of the available males */
               ch = floor(nmale*randunif());
               m  = MALES[ch];
            }
        }
        OFF[i][0] = Ind;                  /* Label for the new individual born */
        Ind++;                            /* Add a new one to Ind */
        OFF[i][1] = floor(2*randunif());  /* Offspring randomly m or f */
        OFF[i][2] = ID[h][2];             /* Mother's nest coordinate x */
        OFF[i][3] = ID[h][3];             /* Mother's nest coordinate y */
        OFF[i][4] = 0;                    /* Starts at age zero */ 
        OFF[i][5] = ID[h][0];             /* Offspring Mother Number */
        OFF[i][6] = ID[m][0];             /* Offspring Father Number */
        MomP = 0; 
        DadP = 0; /* Below finds the matrix pos of Mom, Dad */
        while(OFF[i][5]!=Rmof[MomP][0]){
            MomP++; /* Keep increasing until MomP finds the Mum's position */
        }
        while(OFF[i][6]!=Rmof[DadP][0]){
            DadP++; /* Keep increasing until MomP finds the Mum's position */
        } /* Now have both Mum and Dad position in the kinship matrix `Rmof' */
        OFF[i][7] = Rmof[MomP][DadP+1]; /* Offspring inbreeding coefficient */
        OFF[i][8] = -1;  /* Offspring Mate selection -- no mate selected */
        if(ID[h][8]!=m){
            OFF[i][9] = 1; /* Extra pair offspring */
        }else{
            OFF[i][9] = 0; /* Within pair offspring */
        }
        for(j=0; j<Nloci; j++){ 
            s = floor(randunif()*8); /* Randomly selects a number from 1 to 8 */
            /* There are eight possible ways to pass alleles parent-offspring */
            switch(s){ 
                case 0:
                    OFF[i][10+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][12+(j*4)];
                    break;
                case 1:
                    OFF[i][10+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][13+(j*4)];
                    break;
                case 2:
                    OFF[i][10+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][12+(j*4)];
                    break;
                case 3:
                    OFF[i][10+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[m][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[h][13+(j*4)];
                    break;
                case 4:
                    OFF[i][10+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][12+(j*4)];
                    break;
                case 5:
                    OFF[i][10+(j*4)] = ID[h][10+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][12+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][13+(j*4)];
                case 6:
                    OFF[i][10+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][10+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][12+(j*4)];
                    break;
                case 7:
                    OFF[i][10+(j*4)] = ID[h][11+(j*4)];
                    OFF[i][11+(j*4)] = ID[m][11+(j*4)];
                    OFF[i][12+(j*4)] = ID[h][13+(j*4)];
                    OFF[i][13+(j*4)] = ID[m][13+(j*4)];
                    break;
                default:
                    printf("ERROR IN PARENTS.C\n");
                    break;
            }
        } /* 1 allele of each parent now randomly inherited by offspring at all loci */
    }
    FREE_1ARRAY(MALES); /* Regardless, should be a MALES array to free */
    FREE_1ARRAY(Mse);   /* Regardless, should be a Mse array to free */
    /* ==========================================================*/
    /* Given some conditions, Offspring alleles can mutate ======*/
    /* ==========================================================*/        

    for(j=0; j<l; j++){ /* For all `l' new offspring in OFF */
        for(h=0; h<Nloci; h++){ /* For each offsprings' loci */
            if(randunif() < mu){ /* Copy 1 of loci mutates? */
                a = randnorm(mumu,musd); /* New mutation */
                OFF[j][10+(h*4)] += a; /* Replace old w/new */
            } /* Below does the same as above for 2nd allele */
            if(randunif() < mu){ /* Copy 2 of loci */
                a = randnorm(mumu,musd); /* New mutation */
                OFF[j][11+(h*4)] += a;
            }
        }
    } /* All mutations in the new offspring have now taken place */

    /*===========================================================================*/
    /* PRINTS OFF PEDIGREE BELOW     ============================================*/
    /*===========================================================================*/

    if(prP == 10000){ /* Don't ever want to print this -- unless want FULL pedigree */
        Pedigree = fopen("Pedigree.txt","a+");
        for(i=0; i<l; i++){
            pscore = 0.0;
            nscore = 0.0;
            mscore = 0.0;
            for(j=0; j<Active; j++){
                pscore += OFF[i][((4*j)+10)];
                pscore += OFF[i][((4*j)+11)];
            }
            for(j=0; j<Neutral; j++){
                nscore += OFF[i][((4*j)+Neutstart)];
                nscore += OFF[i][((4*j)+Neutstart+1)];
            }
            for(j=0; j<load; j++){
                mscore += OFF[i][((4*j)+loadstart)];
                mscore += OFF[i][((4*j)+loadstart+1)];
            }


            fprintf(Pedigree,"%d\t%f\t%f\t%d\t%f\t",pid,Beta1,Pcost,gener,OFF[i][0]);
            fprintf(Pedigree,"%f\t%f\t%f\t",OFF[i][1],OFF[i][5],OFF[i][6]);
            fprintf(Pedigree,"%f\t%f\t",OFF[i][7],pscore);
            fprintf(Pedigree,"%f\t%f\t%f\n",nscore,mscore,OFF[i][9]);
        }
        fclose(Pedigree);
    }
}








