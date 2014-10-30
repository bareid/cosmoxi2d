#include "cosmoxi2d.h"
#include "misc.h" //this is where legendrep() lives.
#include "pktosig.h"
#include "pktovin.h"
#include "pktoxi.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

int streamfree();
int dostream(double smin, double smax, double deltas, double isotropicdisp1d, double fvel, int ellmaxin, double alphaperp, double alphapar);
double xistreameval(double sval, int ell);
int printxis(char *outfname, double bias, double fvel, double iso1d, double alphaperp, double alphapar);
