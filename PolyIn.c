/***********************************************************************************/
/* Program: PolyIn.c
   By: Brad Duthie                                             
   Description: Calls Inbreed.c multiple times to simulate different inbreeding
        scenarios.
   Compile: To compile with makefile, type `make' within directory, then hit ENTER
   Run:     To run, type `./PolyIn' on command line, then hit ENTER   */
/***********************************************************************************/
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Inbreed.h"
#include "array.h"

int main(void){

    int    mc, M, Imm, Clu, rep, i, j, xlen, Pedi, prP, conSt, snap, msel, condk;
    int    Active, Neutral, load, gen, muSt, Kind, EpRestr, WpRestr, PreSel, PostSel;
    int    PreserC, wpVsep;
    double Beta1, alpha, mu, *RES, ImmSD, poadj, epadj, wpadj, Scost, Pcost, Ecost;
    double poe, wpe, epe;

    /* =========== VARIABLES BETWEEN THE Xs BELOW ADJUST MODEL PARAMETERS ================*/
    /*  XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX   */
    /* ===================================================================================*/
    /* Model parameter values                                                             */
    /* ===================================================================================*/
    mc      = 100;    /* Number of females a male can mate with (as a social mate)        */
    M       = 0;      /* Maximum age (0 = non-overlapping gen; do not make >0 [ped error])*/
    Imm     = 5;      /* Number of immigrants per generation                              */
    Clu     = 8;      /* Clutch size for each female                                      */
    alpha   = 1.0;    /* Strength of inbreeding avoidance/preference                      */
    Beta1   = 0.0000; /* Selection coefficient on deleterious recessives                  */
    Scost   = 0.0000; /* Cost of having tendancy to non-randomly select social mates      */
    Pcost   = 0.0000; /* Cost of having tendancy to engage in polyandry                   */
    Ecost   = 0.0000; /* Cost of having tendancy to non-randomly select extra-pair mates  */
    gen     = 5000;   /* Number of generations per replicate (default 5000)               */
    muSt    = 0;      /* Generation at which mutations may start                          */
    mu      = 0.001;  /* Mutation rate of any given allele                                */
    Kind    = 1;      /* Kin recognition (1 = recognise all; 0 = recognise only siblings) */
    xlen    = 10;     /* Spatial x and y dimension -- carrying capacity = xlen*xlen       */
    ImmSD   = 1.0;    /* SD around the mean for immigrant allelic values                  */
    conSt   = 1;      /* Constrain WP & EP avoidance/preference to be equal? (0:no, 1:yes)*/
    wpVsep  = 0;      /* 0: wpe is soc, epe is poly mates | 1: wpe is pre, epe is postcop */ 
    condk   = 0;      /* Will females reject EPMs if none better than soc? (1:yes)        */
    EpRestr = 0;      /* Restrict additional mate acces? (0: no, >0: Restriction number   */
    WpRestr = 0;      /* Restrict mate access for WP mates? (0: no, >0: Restriction #     */
    PreSel  = 1;      /* Allow pre-copulatory mate selection for EPMs (0: no, 1:yes)      */
    PostSel = 0;      /* Allow post-copulatory mate selection for EPMs (0: no, 1: yes)    */
    PreserC = 1;      /* Preserve among-individual trait correlations in immigrants       */
    /* ===================================================================================*/
    /* Genome attributes of individuals                                                   */
    /* ===================================================================================*/
    Active     = 10;   /* Number of social pairing alleles (wp preference or avoidance)   */
    Neutral    = 10;   /* Number of polyandry alleles (affect fem taking epms)            */
    load       = 10;   /* Number of extra-paring alleles (ep preference or avoidance)     */
    wpe        = 1.0;  /* Do wp alleles have any effect? (1: yes [default], 0: no)        */
    poe        = 1.0;  /* Do polyandry alleles have any effect? (1: yes [default], 0: no) */
    epe        = 1.0;  /* Do extra-pair alleles have any effect? (1: yes [default], 0: no)*/
    wpadj      =  0.0; /* External adjustment to WP-Pref parameter (default = 0)          */
    poadj      =  0.0; /* External adjustment to Polyandry parameter (default = 0)        */
    epadj      =  0.0; /* External adjustment to EP-Pref parameter (default = 0)          */
    /*                                                                                    */
    /* ===================================================================================*/
    /* Simulation details                                                                 */
    /* ===================================================================================*/
    rep        = 1;     /* Simulations run                                                */
    Pedi       = 0;     /* Print last rep's pedigree? (0:no, 1:yes) WARNING: 200Mb file   */
    snap       = 1;     /* Last 2 gens pedigree for all reps printed? (0:no, 1:yes)       */
    msel       = 1;     /* Last 2 gens print mate selection of females? (0:no, 1:yes)     */
    /*                                                                                    */
    /* ===================================================================================*/
    /*  XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX   */
    /* ===================================================================================*/
    srand(time(NULL) ^ (getpid()<<16)); /* Use time to generate random seed */

    i    = 0; /* Loops through different replicate simulations */
    prP  = 0; /* Indicator for actually printing the pedigree */
    while(i < rep){
        if(i == rep-1 && Pedi == 1){ /* If it's the last rep, and should print pedigree */
            prP = 1;                 /* Switch this variable so pedigree will print */
        }
        MAKE_1ARRAY(RES,11); /* The RES array holds summary statistics */
        for(j=0; j<11; j++){ /* Need to refresh RES with zeros for each simulation */
            RES[j] = 0.0;
        }
        /* The function below is the main simulation function from Inbreed.c */
        Inbreed(mc,M,Imm,Clu,RES,Beta1,i,Active,Neutral,load,alpha,gen,muSt,mu,Kind,xlen,
            prP,ImmSD,wpe,poe,epe,poadj,epadj,wpadj,Scost,Pcost,Ecost,conSt,snap,PostSel,
            msel,EpRestr,WpRestr,condk,PreSel,wpVsep,PreserC);

        FREE_1ARRAY(RES); /* Free the RES array after printing is finished */
        /* Below increases the selection coefficients for the next loop */
        i++;
    }
    return 0;
}
