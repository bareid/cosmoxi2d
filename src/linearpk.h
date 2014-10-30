/***********

Note: units converted to k in Mpc^-1, P(k) in Mpc^3 in getpkspline !!!

***********/

#include "cosmoxi2d.h"
#include "misc.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>  //need this so that gsl_pow will work!!

double pk(double kval);
int getpkspline(char *pkfname, double pkhval, double pkampscale);
