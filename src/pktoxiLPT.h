#include "cosmoxi2d.h"
#include "misc.h"
#include "ptraw.h"
#include "besselroots.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

int getxiLPTmoments(double intkmax, double smin, double smax, double deltas, int ellmax);
int printxiLPTmoments(char *outfname);
double xiLPTmoments(double sval, int ell);
double xiLPTred(double sval, double mu);
int pktoxiLPTfree();
