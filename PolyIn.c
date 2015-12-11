/***********************************************************************************/
/* Program: PolyIn.c
   By: Brad Duthie                                             
   Description: Calls Inbreed.c multiple times to simulate different inbreeding
        scenarios.
   Compile: To compile with makefile, type `make' within directory, then hit ENTER
   Run:     To run, type `./InMult' on command line, then hit ENTER   */
/***********************************************************************************/

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Inbreed.h"
#include "array.h"

int main(void){

    int    mc, M, Imm, Clu, rep, i, j, k, xlen, Pedi, prP, conSt, snap, EpAvail, msel;
    int    Active, Neutral, load, gen, muSt, Kind, mNalleles, sdNalleles, EpRestr, WpRestr;
    int    condk;
    double Beta1, alpha, mu, *RES, ImmSD, poadj, epadj, wpadj, Scost, Pcost, Ecost;
    double poe, wpe, epe;

    /* =========== VARIABLES BETWEEN THE Xs BELOW ADJUST MODEL PARAMETERS ================*/
    /*  XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX   */
    /* ===================================================================================*/
    /* Model parameter values                                                             */
    /* ===================================================================================*/
    mc      = 100;    /* Number of females a male can mate with                           */
    M       = 0;      /* The age at which all individuals must die (0 = non-over)         */
    Imm     = 15;      /* Number of immigrants per generation                              */
    Clu     = 24;      /* Clutch size for each female                                      */
    alpha   = 1.0;    /* Strength of inbreeding avoidance/preference                      */
    Beta1   = 5.0000; /* Selection coefficient on deleterious recessives                  */
    Scost   = 0.0000; /* Cost of having tendancy to non-randomly select social mates      */
    Pcost   = 0.0000; /* Cost of having tendancy to engage in polyandry                   */
    Ecost   = 0.0000; /* Cost of having tendancy to non-randomly select extra-pair mates  */
    gen     = 200;   /* Number of generations per replicate (default 5000)               */
    muSt    = 0;      /* Generation at which mutations may start                          */
    mu      = 0.001;  /* Mutation rate of any given allele                                */
    Kind    = 1;      /* Kin recognition (1 = recognise all; 0 = recognise only siblings) */
    xlen    = 10;     /* Spatial x and y dimension -- carrying capacity = xlen*xlen       */
    ImmSD   = 1.0;    /* SD around the mean for immigrant allelic values                  */
    conSt   = 1;      /* Constrain WP & EP avoidance/preference to be equal? (0:no, 1:yes)*/
    condk   = 0;      /* Is gPy scaled by the kinship of the initial mate (1:yes)         */
    EpRestr = 0;      /* Restrict mate access when EpAvail<=0? (0: no, >0: Restriction #  */
    WpRestr = 0;      /* Restrict mate access for WP mates? (0: no, >0: Restriction #     */
    EpAvail = 0;      /* EP rules: <=0: EP subset selected with EP alleles, random within */
    /*                              >0: EP subset random, EP alleles used within subset   */
    /* ===================================================================================*/
    /* Genome attributes of individuals                                                   */
    /* ===================================================================================*/
    Active     = 10;   /* Number of social pairing alleles (wp preference or avoidance)   */
    Neutral    = 10;   /* Number of polyandry alleles (affect fem taking epms)            */
    load       = 10;   /* Number of extra-paring alleles (ep preference or avoidance)     */
    mNalleles  = 2;    /* Set mean number of alleles per locus (must be >1)               */
    sdNalleles = 0;    /* Set SD of allele number per locus                               */
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
    k    = 0; /* Use for the switch in Beta */
    prP  = 0; /* Indicator for actually printing the pedigree */
    while(i < rep){
        if(i == rep-1 && Pedi == 1){ /* If it's the last rep, and should print pedigree */
            prP = 1;                 /* Switch this variable so pedigree will print */
        }

        MAKE_1ARRAY(RES,11); /* The RES array holds summary statistics */
        for(j=0; j<11; j++){ /* Need to refresh RES with zeros for each simulation */
            RES[j] = 0.0;
        }
        /* XXX XXX <<<< MAKE VECTOR 1:1000ish TO ORDER SNAP PED OUTPUTS LATER >>>>> XXX */
        /* The function below is the main simulation function from Inbreed.c */
        Inbreed(mc,M,Imm,Clu,RES,Beta1,i,Active,Neutral,load,alpha,gen,muSt,mu,Kind,xlen,
            mNalleles,sdNalleles,prP,ImmSD,wpe,poe,epe,poadj,epadj,wpadj,Scost,Pcost,Ecost,
            conSt,snap,EpAvail,msel,EpRestr,WpRestr,condk);
 
        FREE_1ARRAY(RES); /* Free the RES array after printing is finished */
        /* Below increases the selection coefficients for the next loop */
        i++;
        switch(k){
            case 0:
                Beta1 = 0.2;
                k++;
                break;
            case 1:
                Beta1 = 1;
                k++;
                break;
            case 2:
                Beta1 = 2;
                k++;
                break;
            case 3:
                Beta1 = 5;
                k++;
                break;
            default:
                Beta1 = 0;
                Pcost += 0.01;
                k = 0;
                break;
        }
    }
    return 0;
}

