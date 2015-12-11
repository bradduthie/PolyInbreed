#include <stdio.h>
#include <math.h>
#include "array.h"


void RegrTr(double **REGR, double **OFF, int prevOff);

void inbrdep(double **ID, double **REGR, int prevOff, int Liv, int M);

void printMeanKinship(double **Rmof, double **ID, int Liv,
    int Active, int Neutral, int Neutstart, int load, int loadstart, int pid,
    double Beta1, double Pcost, int WpRestr, int EpRestr);
