/***********************************************************************************/
/* Program: Inbreed.c
   By: Brad Duthie                                             
   Description: IBM to simulate inbreeding tolerance/avoidance
   Compile: gcc Inbreed.c -ansi -Wall -pedantic    */
/***********************************************************************************/

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "array.h"
#include "randnormINT.h"
#include "randnorm.h"
#include "landscape.h"
#include "arraysort2D.h"
#include "vectorrank.h"
#include "Rbuild.h"
#include "mateselect.h"
#include "babies.h"
#include "parents.h"
#include "mortality.h"
#include "retain.h"
#include "deadbottom.h"
#include "sumstats.h"
#include "inbrdep.h"
#include "SampleNoReplace.h"

#define CORES 1 /* For running in parallel, set cores > 1 */
#define length(x) (sizeof (x) / sizeof *(x))

void Inbreed(int mc, int M, int Imm, int Clu, double *RES, double Beta1, int rep,
    int Active, int Neutral, int load, double alpha, int gen, int muSt, double mu, 
    int Kind, int xlen, int prP, double ImmSD, int wpe, int poe, int epe, 
    double poadj, double epadj, double wpadj, double Scost, double Pcost, 
    double Ecost, int conSt, int snap, int PostSel, int msel, int EpRestr, 
    int WpRestr, int condk, int PreSel, int wpVsep, int PreserC){

    /* =========================================================================*/
    /* =========================================================================*/
    /* DEFINE VARIABLES TO BE USED: */
    /* =========================================================================*/
    /* =========================================================================*/

    int i, j, k, l, h, g; 
    int Ind, suc, Liv, OLiv, cols, prevOff, removed, pid, pidadd, pidcmb;
    int Nloci;         /* The number of loci in simulated genome */
    int TotAlleles;    /* Keep track of total alleles in sim */
    int StInds;        /* The initial number of females and males */
    int *O;            /* Keeps track of offspring per nest in a generation */
    int loadstart;     /* where inbreeding load loci start on table */
    int Neutstart;     /* Where do the neutral's start on the table */
    int *shuffImm;     /* Shuffle vector for immigrants */

    double a;          /* Scalar place for alleles that can mutate */
    double a1[2];      /* Vector used to transition making alleles */
    double *alleles;   /* A vector used in generating alleles */
    double BETA;       /* Any habitat autocorrelation -- none is 0; more is < 0 */
    double *LS;        /* Used to build a landscape */
    double **SCAPE;    /* This is the landscape itself */
    double **GenArch;  /* 2D array to hold genetic architecture */
    double **ID;       /* 2D array to temporarily hold individuals and traits */
    double **Rank;     /* lists quality rank of each x y landscape coordinate */ 
    double **Rmat;     /* Coefficient of relatedness matrix for parents*/
    double **Rmof;     /* Coefficient of relatedness matrix for offspring */
    double **OFF;      /* 2D array to temporarily hold offspring */
    double **PIV;      /* An array to combine ID & OFF and pivot to one ID */
    double **REGR;     /* An array to hold values for regression calcs */
    double musd;       /* Standard deviation in mutational effects */
    double mumu;       /* Mean mutational effect */
    double *rImm;      /* Random value for an immigrant */

    char SnapPed[20];  /* Snapshot pedigree -- last generation recorded */
    char evolres[20];  /* Evolution over multiple generations */

    FILE *fptr;        /* Snap shot pedigree file (see above) */
    FILE *evol;        /* Evolution over multiple generations (see above) */

    /* =========================================================================*/
    /* =========================================================================*/
    /* INITIALISATION:                                                          */
    /* =========================================================================*/
    /* =========================================================================*/

	pid    = getpid(); /* Print out the unique number (for parallel) */

    pidadd = floor(randunif()*10000); /* A bit of insurance for printing file */

    i = 1; /* This just makes a combined number for file output on the cluster */
    while (i <= pidadd){ /* So if running a bunch in parallel, it will be even */
        i *= 10; /* less likely that two file names will print out the same */
    } /* (this would require the process and random number to be identical) */
    pidcmb = pid*i + pidadd; 

    Nloci  = Active + Neutral + load; /* Set the number of loci */

    StInds = 100;   /* Starting number of females and males */
    Ind    = 0;     /* Running count of the number of individuals */
    Liv    = 0;     /* Running count of number of live individuals */

    loadstart = 10 + (4*Active) + (4*Neutral); /* First load allele starts */
    Neutstart = 10 + (4*Active);               /* First neutral allele starts */

    BETA = 0.0;            /* Spatial autocorrelation on the landscape if > 0 */

    cols    = (10+(4*Nloci)); /* Number of cols in the big matrix */
    prevOff = 0; /* No offspring in the previous generation at the start */
    MAKE_2ARRAY(REGR,StInds,8); /* An array to be freed, then used later */

    musd    = sqrt((1.0 / 2 * Active));
    mumu    = 0.0;

    /* =========================================================================*/
    /* Generate genetic architecture ===========================================*/
    /* =========================================================================*/
    TotAlleles = 0;              /* Start counting with zero total alleles */
    MAKE_1ARRAY(alleles, Nloci); /* Make a vector for alleles per loci */
    a = 0;                       /* Use 'a' here to check allele count is > 0. */
    for(i=0; i<Nloci; i++){      /* Loop assigns allele number to loci (def = 2) */
        a = randnormINT(2, 0);
        if(a < 1){ /* Need at least one allele per loci */
            a = 1;
        }
        alleles[i] = a;
        TotAlleles = TotAlleles + a; /* Add a to the total number of alleles */
    } /* Now `alleles' has allele number for each loci, total alleles also counted */
    MAKE_2ARRAY(GenArch,TotAlleles,3); /* Keeps genetic architecture in an array*/
    a = 0;    /* Keep track of allele number in making table */
    l = 0;  /* Keep track of locus number in making table */
    /* Below makes the table showing loci (row 0), allele (1), and allele value (2) */
    for(i=0; i<TotAlleles; i++){
        GenArch[i][0] = l; /* Row zero is denotes the locus */
        GenArch[i][1] = a; /* Row one denotes the allele */
        GenArch[i][2] = randnorm(0,1); /* Row two denotes some allele value */ 
        if((alleles[l]-1)==a){ /*Need -1 bc count one less than label (a) */ 
            l++;
            a = 0;
        }else{
            a++;    
        }
    } /* End result gives a table to look up any allele at any loci */
    /* Note, alllele values are not used at the moment -- and should not be */
    /* Especially for the first position alleles (strategy) */

    /* =========================================================================*/
    /* Generate females and males ============================================= */
    /* =========================================================================*/
    a1[0] = 0; /* Will need an alleles dummy variable here within below loops */
    a1[1] = 0;
    MAKE_2ARRAY(ID,StInds*2,cols); /* This is the central array holding individuals */
    for(i=0; i<(StInds*2); i++){   /* The loop adds starting individual traits */
        ID[i][0] = Ind; /* Individual Number */
        if(i<StInds){   /* Makes StInds females and StInds males in order */
            ID[i][1] = 0;  /* Female */
        }else{
            ID[i][1] = 1;  /* Male */
        }
        ID[i][2] = floor(randunif()*xlen); /* rand x location */
        ID[i][3] = floor(randunif()*xlen); /* rand y location */
        ID[i][4] = 0;                      /* Age -- -1 indicates dead */
        ID[i][5] = -1;                     /* Mother (-1s indicate immigrant) */
        ID[i][6] = -1;                     /* Father (-1s indicate immigrant) */
        ID[i][7] = 0;                      /* Coefficient of inbreeding (F) */
        ID[i][8] = -1;                     /* Mate selection (-1 no mate yet) */
        ID[i][9] = 0;                      /* Within (0) or extra-pair (1) ? */
        for(j=0; j<Nloci; j++){ /* This loop builds starting individual genomes */
            /* Note below, all 0 alleles are made at the start */
            a1[0] = floor(randunif()*(alleles[j]-1));
            a1[1] = floor(randunif()*(alleles[j]-1));
            ID[i][((4*j)+10)] = a1[0]; /* Adds first allele at locus */
            ID[i][((4*j)+11)] = a1[1]; /* Adds second allele at locus */
            suc = 0; /* Triggers break below when allele values are added */
            for(k=0; k<TotAlleles; k++){ /*Adds norm distr. vals */
                if(j==GenArch[k][0] && a1[0]==GenArch[k][1]){
                    ID[i][((4*j)+12)] = GenArch[k][2];
                    suc++; 
                }
                if(j==GenArch[k][0] && a1[1]==GenArch[k][1]){
                    ID[i][((4*j)+13)] = GenArch[k][2];
                    suc++;
                }
                if(suc==2){ /* Avoids looping needlessly */
                    break; /* And breaks when both alleles are found */
                } /* The individual is now complete */
            }
        }
        Ind++; /* Move on to the next individual number */
    }
    Liv = Ind; /*Keep track of total living individuals */

    /* =========================================================================*/
    /* Generate nest sites/landscape ========================================== */
    /* =========================================================================*/
    /* All of this builds an explicit landscape, if there is reason to use one */
    /* Recall below is same: LS = (double*)malloc(sizeof(double)*(xlen*xlen)); */
    MAKE_1ARRAY(LS,(xlen*xlen));
    MAKE_2ARRAY(SCAPE,xlen,xlen);
    MAKE_2ARRAY(Rank,xlen*xlen,3); /* Keeps genetic architecture */

    /* The landscape LS will be created below */
    landscape(LS,BETA,xlen); /* Autocorrelation dependent upon BETA */
    /* Below, the each cell value from above is fed into SCAPE */
    for(i=0; i<xlen; i++){
        for(j=0; j<xlen; j++){
            SCAPE[i][j] = LS[(xlen*i)+j];
        }
    } /* Now each landscape cell has a value, which can be used if needed */

    /* =========================================================================*/
    /* Generate a dummy r matrix to start ===================================== */
    /* =========================================================================*/
    
    MAKE_2ARRAY(Rmat, Liv, Liv+1); /* All live, first col head */

    for(i=0; i<Liv; i++){ /* Adds the ind number to Rmat col 0 */
        Rmat[i][0] = -1;
    } /* Now Rmat is not square, but additional col for header */

    for(i=0; i<Liv; i++){ 
        for(j=0; j<Liv; j++){
            if(j!=i){ /* 0s for all off-diagonal */
                Rmat[j][i+1] = 0;
            }else{ /* Diagonial is always r = 1 to start*/
                Rmat[j][i+1] = 1;
            }
        }
    }

    OLiv = Liv; /* Needed in Rbuild when Rmat and Rmof matrix sizes differ. */
    /* OLiv will be the number of individuals in the old matrix, Liv in the new */

    /* =========================================================================*/
    /* =========================================================================*/
    /* BEGIN SIMULATION RUN                                                     */
    /* =========================================================================*/
    /* =========================================================================*/

    i = 0;
    while(i < gen){
        /* ==========================================================*/
        /* See if population size is below extinction threshold      */
        /* ==========================================================*/    
        g = 0; /* Will use to count the number of living females */
        k = 0; /* Will use to count the number of living males   */
        for(j=0; j<Liv; j++){
            if(ID[j][4] <= M && ID[j][1] == 0){
                g++; /* Adds to `g' if alive, or recently dead */
            } /* And if female */
            if(ID[j][4] <= M && ID[j][1] == 1){
                k++; /* Adds to `k' if alive, or recently dead */
            } /* And if male */
        } /* If male or female number is less than immigrants, extinct */
        if(g <= Imm || k <= Imm){ 
            RES[0] = -1;
            RES[1] = -1;
            break; /* This busts out of the while loop -- ends simulation */
        }
        /* ==========================================================*/
        /* Randomise & rank individuals                              */
        /* ==========================================================*/    
        for(j=0; j<Liv; j++){
            ID[j][8]  = -1; /* Ensures no mate selected yet */
            ID[j][loadstart+2] = randunif(); /* Random number in allele 1 value place */
            if(ID[j][1]==0){ /* Random number used to rank individuals */
                /* Shift down 2 if female (need at top after sort) */
                ID[j][loadstart+2] = ID[j][loadstart+2] - 2; 
                /* Given tendancy to non-randomly mate, chance of dying */
                if(wpe == 1 && Scost > 0){ /* Within-pair strategy alleles effect? */
                    g = 0;
                    for(h=0; h<Active; h++){
                        g += ID[j][((4*h)+10)];
                        g += ID[j][((4*h)+11)];
                    }  /* no change if poe=0 (polyandry alleles inactive) */
                    g *= Scost; /* How great tendancy for polyandry? */
                    if(g < 0){
                        g *= -1.0;
                    }
                    if(g > randunif()){ /* Kill from cost ? */
                        ID[j][4] = -1;
                        ID[j][loadstart+2] = ID[j][loadstart+2] + 90;
                    }
                }
                /* Given tendancy for polyandry, chance of dying */
                if(poe == 1 && Pcost > 0){ /* If polyandry alleles have an effect */
                    g = 0;
                    for(h=0; h<Neutral; h++){
                        g += ID[j][((4*h)+Neutstart)];
                        g += ID[j][((4*h)+Neutstart+1)];
                    }  /* no change if poe=0 (polyandry alleles inactive) */
                    g *= Pcost; /* How great tendancy for polyandry? */
                    /* Note: If g < 0, no cost because no tendancy for polyandry */
                    if(g > randunif()){ /* Kill from cost ? */
                        ID[j][4] = -1;
                        ID[j][loadstart+2] = ID[j][loadstart+2] + 90;
                    }
                } /* be effectively out of the pool for finding a mate */
                if(epe == 1 && Ecost > 0){ /* EP strategy alleles effect? */
                    g = 0;
                    for(h=0; h<load; h++){
                        g += ID[j][((4*h)+loadstart)];
                        g += ID[j][((4*h)+loadstart+1)];
                    }  /* no change if poe=0 (polyandry alleles inactive) */
                    g *= Ecost; /* How great tendancy for polyandry? */
                    if(g < 0){
                        g *= -1.0;
                    }
                    if(g > randunif()){ /* Kill from cost ? */
                        ID[j][4] = -1;
                        ID[j][loadstart+2] = ID[j][loadstart+2] + 90;
                    }
                }
            }
            if(ID[j][4]<0 || ID[j][4]>M){ /* All dead way down, pre-sort */
                ID[j][loadstart+2] = ID[j][loadstart+2] + 10.0;
            } /* The +10 makes dead individuals on the bottom when array sorted */
        } /* Females now lowest values in col 12, deads all have highest -- can sort */
        arraysort2D(ID, Liv, cols, loadstart+2, 1);  /* Sorts by rand unif in col ldst+2 */ 

        /* ==========================================================*/
        /* Females compete for nest sites */
        /* ==========================================================*/    

		/* First determine sites would have been held by individuals lost to the cost */
        removed = 0; /* This will figure out how many removed due to a direct cost... */
        for(j=0; j<Liv; j++){ /* that WOULD have otherwise made it into the mating pool */
            if(ID[j][loadstart+2] > 50){ /* A rand in col 12 >50 happens iff cost realised */
                k = 0; /* So set the count to zero */
                while((ID[j][loadstart+2]-100)>ID[k][loadstart+2] 
                       && k<=xlen*xlen-removed && ID[k][1]==0){
                    k++; /* If the 100 (90 for cost, 10 for dead) subtracted ... */
                } /* put it in the top xlen*xlen (corrected for already removed)... */
                if(k <= xlen*xlen){ /* then record this as has been removed  */
                    removed++; /* Another one bites the dust */
                } /* This effectivley makes it possible to make it as though inds...*/
            } /* Were removed after applying carrying capacity, instead of before */
        } /* This second for j loop avoids two computation-heavy array sorts */

        /* Below ranks the nest sites by value in an easy lookup */
        vectorrank(LS, Rank, xlen); /* Random unless landscape explicit */

        j = 0; /* This assigns females to nest */
        k = 0; /* If not enough nests, female pos is -1 */
        while(j < Liv && k < (xlen*xlen - removed)){ /* As long as sites exist: */
            /* If an individual is alive and female */
            if(ID[j][4] >= 0 && ID[j][4] <= M && ID[j][1]==0){
                ID[j][2] = Rank[k][0]; /* Give it an x location */
                ID[j][3] = Rank[k][1]; /* Give it an y location */
            }else{                     /* If it is not alive */
                ID[j][2] = -1;         /* Don't give it any position */
                ID[j][3] = -1;
            }
            j++; /* Keeps going until the array ID rows are depleted */
            k++; /* Or until the number of sites available are depleted */
        }
        if(Liv > (xlen*xlen - removed)){ /* If more inds in ID array than territories */
            for(j=k; j<Liv; j++){  /* Start where the territories run out */
                ID[j][2] = -1; /* Then remove the locations from any inds */
                ID[j][3] = -1;
            }
        } /* All able choosing sex individuals should now have territories */

        /* ==========================================================*/
        /* Cull non-choosing sex  -- by forcing to find a site too   */
        /* ==========================================================*/    
        k = 0; /* Note that this forces a 1:1 sex ratio if > carrying capacity */
        while(ID[k][1] == 0 && k < Liv){
            k++; /* k is where the non-choosing sex starts */
        }
        j = 0; /* Loop below goes through all non-choosing sex, assigns territory */
        while(ID[k][1] != 0 && j < (xlen*xlen) && j < Liv-k){ 
            ID[k][2] = Rank[j][0]; /* If available (above while condition) */
            ID[k][3] = Rank[j][1];
            j++;
            k++;
        } /* Choosing sex will later only mate with ind on a site */

        if(i > 0){ /* If we're not looking at the first generation */
            inbrdep(ID, REGR, prevOff, Liv, M); /* Add ending survival ? to REGR */
        }
        /* ===============================================================*/
        /* Generate array of kinship coefficients from ID parenthood data */
        /* ===============================================================*/

        MAKE_2ARRAY(Rmof, Liv, Liv+1); /* All live, first col head */

        for(j=0; j<Liv; j++){ /* Adds the ind number to Rmat col 0 */
            Rmof[j][0] = ID[j][0]; 
        } /* Rmof is not square, but additional col for header */

        /* Rbuild makes a new R matrix for living  based on previous R */
        Rbuild(ID, Rmat, Rmof, Liv, OLiv);

        if(snap == 2 && i == (gen - 1)){
            printMeanKinship(Rmof, ID, Liv, Active, Neutral, Neutstart, load, loadstart,
                pidcmb, Beta1, Pcost, WpRestr, EpRestr);
        } /* This prints out all kinship coefficients in the last generation */

        /* ==========================================================*/
        /* Choosing sex competes for and chooses pair bonded mates   */
        /* ==========================================================*/        

        mateselect(ID, Rmof, Liv, Nloci, M, mc, Active, alpha, Kind, wpe, wpadj, WpRestr);

        /* ==========================================================*/
        /* Determine the number of offspring in each nest            */
        /* ==========================================================*/        

        k = 0; /*Below counts the number of the females */ 
        for(j=0; j<Liv; j++){ /* Must be female, have mate, be alive, nest */
            if(ID[j][1]==0 && ID[j][8] > -1 && ID[j][4] > -1 
                && ID[j][2]>-1 && ID[j][4] <= M){
                k++;
            }else{
                break; /* Ordered, so stop if not female with nest */
            } /* k is now the number of choosing sex individuals breeding */
        } 
        for(j=0; j<Liv; j++){ /* Returns males from -1 values */
            if(ID[j][1]<0){ /* Were negative from mateselect */
                ID[j][1] = 1; /* This brings them back */
            }
        } 

        MAKE_1ARRAY(O, k); /* To measure each nest's offspring number */

        l = 0;
        for(j=0; j<k; j++){ /* Babies adds offspring number for each chooser */
            O[j] = babies(ID, Clu, j) + l; /* Each element sums prev. */
            l    = O[j]; /* l is the total new offspring produced */
        } /* Now have total offspring for each choooser (j) in vector O */

        MAKE_2ARRAY(OFF,l,cols); /* Offspring array (OFF) with l rows (number new): */

        /* ==========================================================*/
        /* Determine who the parents are of each offspring           */
        /* ==========================================================*/        

        parents(ID,OFF,Rmof,O,Nloci,Ind,k,l,Clu,Liv,i,Beta1,alpha,RES,rep,loadstart,
            load,Active,prP,Neutstart,Neutral,Kind,M,poe,epe,poadj,epadj,mu,mumu,
            musd,conSt,PostSel,pidcmb,msel,gen,EpRestr,Pcost,condk,PreSel,wpVsep);
        Ind += l; /*Needed because the Ind++ in parents wont link here */

        /* ==========================================================*/
        /* If last gen, prints last 2 gen pedigree to file           */
        /* ==========================================================*/        
        if(snap == 1 && i == (gen - 1)){
            sprintf(SnapPed,"ped%d.txt",pidcmb);
            fptr = fopen(SnapPed,"a+");
            for(j=0; j<Liv; j++){ /* Printing relevant pedigree information */
                if(ID[j][4]==0){ /* Only concerned with living individuals */
                    fprintf(fptr,"%d\t%f\t%f\t%f\t",pidcmb,ID[j][0],ID[j][1],ID[j][4]+1);
                    fprintf(fptr,"%f\t%f\t%f\t%f\t",ID[j][5],ID[j][6],ID[j][7],ID[j][8]);
                    fprintf(fptr,"%f\t"            ,ID[j][9]);
                    a = 0.0; /* Just to add up the WP values, not ind alleles */
                    for(g=0; g<Active; g++){
                        a += ID[j][((4*g)+10)];
                        a += ID[j][((4*g)+11)];
                    }
                    fprintf(fptr,"%f\t",a);
                    a = 0.0; /* Just to add up the Poly values, not ind alleles */
                    for(g=0; g<Neutral; g++){
                        a += ID[j][((4*g)+Neutstart)];
                        a += ID[j][((4*g)+Neutstart+1)];
                    }
                    fprintf(fptr,"%f\t",a);
                    a = 0.0; /* Just to add up the EP values, not ind alleles */
                    for(g=0; g<load; g++){
                        a += ID[j][((4*g)+loadstart)];
                        a += ID[j][((4*g)+loadstart+1)];
                    }
                    fprintf(fptr,"%f\t%f\t%f\n",a,Beta1,Pcost);
                }
            }
            for(j=0; j<l; j++){ /* Only concerned here with fresh offspring */
                fprintf(fptr,"%d\t%f\t%f\t%f\t",pidcmb,OFF[j][0],OFF[j][1],OFF[j][4]);
                fprintf(fptr,"%f\t%f\t%f\t%f\t",OFF[j][5],OFF[j][6],OFF[j][7],OFF[j][8]);
                fprintf(fptr,"%f\t"            ,OFF[j][9]);
                a = 0.0; /* Just to add up the WP values, not ind alleles */
                for(g=0; g<Active; g++){
                    a += OFF[j][((4*g)+10)];
                    a += OFF[j][((4*g)+11)];
                }
                fprintf(fptr,"%f\t",a);
                a = 0.0; /* Just to add up the Poly values, not ind alleles */
                for(g=0; g<Neutral; g++){
                    a += OFF[j][((4*g)+Neutstart)];
                    a += OFF[j][((4*g)+Neutstart+1)];
                }
                fprintf(fptr,"%f\t",a);
                a = 0.0; /* Just to add up the EP values, not ind alleles */
                for(g=0; g<load; g++){
                    a += OFF[j][((4*g)+loadstart)];
                    a += OFF[j][((4*g)+loadstart+1)];
                }
                fprintf(fptr,"%f\t%f\t%f\n",a,Beta1,Pcost);
            }
            fclose(fptr);
        }

        /* ==========================================================*/
        /* Now some individuals in both ID and OFF die ==============*/
        /* ==========================================================*/    

        /* First ID individuals have aged 1 generation */
        for(j=0; j<Liv; j++){   /* Everyone ages, unless already dead */
            if(ID[j][4] > -1){  /* If they are not dead */
                ID[j][4]++;     /* One year added here */
            }else{              /* If they are dead */
                ID[j][4] = -1;  /* They stay dead */
            }
        } /* All the old individuals from ID now aged/dead */

        /* Now adults (if too old) or juveniles (if too del. recessive) die */
        mortality(ID,OFF,Liv,l,M,Beta1); 

        /* ==========================================================*/
        /* Collect summary statistics                                */
        /* ==========================================================*/        

        sumstats(RES,OFF,M,Nloci,loadstart,load,Neutstart,
            Neutral,Active,l,Imm,xlen,REGR,i,prevOff);
        
        FREE_2ARRAY(REGR);     /* Kills old REGR after sumstats calculations */
        MAKE_2ARRAY(REGR,l,11); /* Some lifetime values for the regression */

        RegrTr(REGR, OFF, l);  /* Transfers relevant info from OFF to REGR */
        prevOff = l;           /* Will need to use this later */

        /* ==========================================================*/
        /* Living, or dead parents with living offspring retained ===*/
        /* ==========================================================*/        

        retain(ID, OFF, Liv, l, M); /* Keeps parents of living to calculate kinships */

        h = 0; /* This keeps track of how many individuals go in a pivot array */
        for(j=0; j<Liv; j++){ /* If an individual is not dead, they get moved  */
            if(ID[j][4] != -1){
                h++; /* Add to h to get another row for the new array */
            }
        }
        for(j=0; j<l; j++){ /* Here now for the offspring array */
            if(OFF[j][4] != -1){ /* We check to see if the offspring is dead */
                h++; /* If not, another row will be needed in the new array */
            }
        }

        h += Imm; /* If there are immigrants to add, also needed in new array */

        MAKE_2ARRAY(PIV,h,cols); /* This matrix pivots ID with all needed Inds */

        h = 0; /* Stick the old live/parent individuals into the pivot */
        for(j=0; j<l; j++){ /* Stick the offspring into the pivot */
            if(OFF[j][4] != -1){
                for(g=0; g<cols; g++){
                    PIV[h][g] = OFF[j][g];
                }                
            h++;
            }
        }
        FREE_2ARRAY(OFF); /* No longer need the offspring matrix */

        for(j=0; j<Liv; j++){ /* Stick the remaining old individuals in the pivot */
            if(ID[j][4] != -1){
                for(g=0; g<cols; g++){
                    PIV[h][g] = ID[j][g];
                }                
            h++;
            }
        } /* Now we have an array PIV with relevant individuals for new gen */

        /* ==========================================================*/
        /* If there are immigrants, we bring them in here ===========*/
        /* ==========================================================*/        

        a1[0] = 0; /* Will need an alleles variables again */
        a1[1] = 0;
        for(j=h; j<(h + Imm); j++){
            PIV[j][0] = Ind; /* Individual Number */ 
            Ind++;           /* Increase Ind to make next individual unique */
            PIV[j][1] = 1;   /* All immigrants male */
            PIV[j][2] = floor(randunif()*xlen); /* rand x location */
            PIV[j][3] = floor(randunif()*xlen); /* rand y location */
            PIV[j][4] = 0;   /* Aged one year on immigration */
            PIV[j][5] = -1;  /* Mother (-1s indicate immigrant) */
            PIV[j][6] = -1;  /* Father (-1s indicate immigrant) */
            PIV[j][7] = 0;   /* Coefficient of inbreeding (F) */
            PIV[j][8] = -1;  /* Mate selection */
            PIV[j][9] = -1;  /* Neigher WP or EP offspring */
            if(PreserC == 0 || i<1){
                for(g=0; g<Nloci; g++){ /* Maintains allele freqs */
                    /* If past the generation of mutation, let */
                    if(g < Active){ /* If it's the allele affecting social mate choice */
                        a1[0] = randnorm(RES[0],ImmSD*RES[3]);
                        a1[1] = randnorm(RES[0],ImmSD*RES[3]);
                    } /* If allele affecting degree of polyandry */
                    if(g >= Active && g < (Active+Neutral)){
                        a1[0] = randnorm(RES[1],ImmSD*RES[4]);
                        a1[1] = randnorm(RES[1],ImmSD*RES[4]);
                    } /* If allele affecting extra pair mate choice */
                    if(g >= (Active+Neutral)){
                        a1[0] = randnorm(RES[2],ImmSD*RES[5]);
                        a1[1] = randnorm(RES[2],ImmSD*RES[5]);
                    }
                    PIV[j][((4*g)+10)] = a1[0]; /* Adds the allele to PIV */
                    PIV[j][((4*g)+11)] = a1[1]; /* Adds the allele to PIV */
                }
            }else{ /* This begins what happens when we preserve correlated traits in */
                if(Neutral != Active){ /* immigrants coming into the population */
                    printf("ERROR: NEUTRAL & ACTIVE *REALLY* SHOULD BE EQUAL IF PreserC");
                } /* What we're doing here is making a lot of random uniforms in advance */
                MAKE_1ARRAY(rImm, cols);
                for(g=0; g<cols; g++){
                    rImm[g] = randnorm(0,ImmSD);
                } /* A shuffling needs to occur to prevent *loci* from being correlated */
                MAKE_1ARRAY(shuffImm, Active); /* The end is that traits are correlated */
                g = 0; /* As they are in the native population, but the shuffle makes */
                while(g < Active){ /* individual loci independent (e.g., L1 for */
                    shuffImm[g] = g;  /* inbreeding and L1 for polyandry ) */
                    g++; /* This has been double-checked to confirm the pattern of */
                } /* output is as desired -- sd within loci, and cor among traits */
                shufflesamp(shuffImm, Active);
                for(g=0; g<Nloci; g++){ /* Maintains allele freqs */
                    /* If past the generation of mutation, let */
                    if(g < Active){ /* If it's the allele affecting social mate choice */
                        a1[0] = rImm[(4*shuffImm[g])+10] * RES[3] + RES[0];
                        a1[1] = rImm[(4*shuffImm[g])+11] * RES[3] + RES[0];
                    } /* Below is icky, but it works (confirmed statistically in R) */
                    if(g >= Active && g < (Active+Neutral)){
                        a1[0]   = sqrt(1-RES[7]*RES[7])*RES[4] * rImm[(4*g)+10] +
                                      RES[7]*RES[4]*rImm[(4*g)-(4*Active)+10] + RES[1];
                        a1[1]   = sqrt(1-RES[7]*RES[7])*RES[4] * rImm[(4*g)+11] +
                                      RES[7]*RES[4]*rImm[(4*g)-(4*Active)+11] + RES[1];
                    } /* If allele affecting extra pair mate choice */
                    if(g >= (Active+Neutral)){
                        a1[0] = rImm[(4*g)+10]  * RES[5] + RES[2];
                        a1[1] = rImm[(4*g)+11]  * RES[5] + RES[2];
                    }
                    PIV[j][((4*g)+10)] = a1[0]; /* Adds the allele to PIV */
                    PIV[j][((4*g)+11)] = a1[1]; /* Adds the allele to PIV */
                }
                FREE_1ARRAY(shuffImm);
                FREE_1ARRAY(rImm);
            } /* Finished what happens if immigrants have correlated traits. */
            Ind++; /* Move on to the next individual number */
        }
 
        /* ==========================================================*/
        /* Remake the matrix ID =====================================*/
        /* ==========================================================*/        

        OLiv = Liv;      /* Keep a hold of the old Liv value, needed later */
        Liv = h + Imm;   /* New Liv value from row number in PIV */
        FREE_2ARRAY(ID); /* We can now get rid of ID */
        MAKE_2ARRAY(ID, Liv, cols); /* But remake it with PIV */
        for(j=0; j<Liv; j++){
            for(g=0; g<cols; g++){
                ID[j][g] = PIV[j][g]; /* Replaces PIV with ID */
            }                
        }        
        FREE_2ARRAY(PIV); /* Then free PIV  */

        /* ==========================================================*/
        /* Now we need to make a new R matrix... ====================*/
        /* ==========================================================*/        

        FREE_2ARRAY(Rmat); /* Get's rid of the old Rmat -- need a new */

        MAKE_2ARRAY(Rmat, OLiv, OLiv+1);

        for(j=0; j<OLiv; j++){ /* Replaces Rmat with Rmof */
            for(k=0; k<(OLiv+1); k++){
                Rmat[j][k] = Rmof[j][k];
            }
        }

        FREE_2ARRAY(Rmof);
        FREE_1ARRAY(O);

        /* ==========================================================*/
        /* Prints generation info in terminal    ====================*/
        /* ==========================================================*/

        sprintf(evolres,"evo%d.txt",pidcmb);
        evol = fopen(evolres,"a+");
        fprintf(evol,"%d\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",pidcmb,
            Beta1,Pcost,i,RES[0],RES[1],RES[2],RES[3],RES[4],RES[5],RES[6],RES[7]);
        fclose(evol);

        /* ==========================================================*/
        /* Move ahead to the next generation                         */
        /* ==========================================================*/        

        i++; /* Now have fresh Rmat and ID arrays to start the next generation */
    }

    FREE_2ARRAY(REGR);
    FREE_2ARRAY(Rmat);
    FREE_2ARRAY(ID);
    FREE_1ARRAY(LS);
    FREE_1ARRAY(alleles);
    FREE_2ARRAY(SCAPE);
    FREE_2ARRAY(GenArch);
    FREE_2ARRAY(Rank);
}


