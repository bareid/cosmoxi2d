#include "cosmoxi2d.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>  //if you DON'T include this, does WACK things but doesn't give an error for gsl_pow_n etc!!! ACK!

#define BMAX 1500

int fillbesselroots();
double getbesselroot(int whichroot, int ell);

int fillbesselrootsxiLPT(int ellmax);
int printbesselrootsxiLPT(int ellmax);
double getbesselrootxiLPT(int whichroot, int ell);

//other oscillatory functions live here.
double mybessel(int whichint, double barg); //used in pktosig.c
double mybesselA13(int whichint, double barg); //oscillatory functions from A13
double getbesselrootA13(int whichroot, int whichbessel);
int fillbesselrootsA13();

