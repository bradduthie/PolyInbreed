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

#define CORES 2
#define length(x) (sizeof (x) / sizeof *(x))

void Inbreed(int mc, int M, int Imm, int Clu, double *RES, double Beta1, int rep,
    int Active, int Neutral, int load, double alpha, int gen, int muSt, double mu, 
    int Kind, int xlen, int prP, double ImmSD, int wpe, int poe, int epe, 
    double poadj, double epadj, double wpadj, double Scost, double Pcost, 
    double Ecost, int conSt, int snap, int PostSel, int msel, int EpRestr, 
    int WpRestr, int condk, int PreSel);
