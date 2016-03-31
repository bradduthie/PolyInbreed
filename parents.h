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
    double mu, double mumu, double musd, int conSt, int PostSel, int pid, int msel,
    int lastgen, int EpRestr, double Pcost, int condk, int PreSel);
