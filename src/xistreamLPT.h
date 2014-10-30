#include "cosmoxi2d.h"
#include "misc.h" //this is where legendrep() lives.
#include "pktoxiLPT.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

int streamLPTfree();
int dostreamxiLPT(double smin, double smax, double deltas, double isotropicdisp1d, int ellmaxin, double alphaperp, double alphapar);
double xistreamLPTeval(double sval, int ell);
int printxisLPT(char *outfname, double bias, double fvel, double iso1d, double alphaperp, double alphapar);

