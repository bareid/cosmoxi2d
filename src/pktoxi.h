#include "cosmoxi2d.h"
#include "linearpk.h"
#include "ptraw.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "besselroots.h"

int getxireal(double intkmax, double rmin, double rmax, double deltar);
int printxireal(char *outfname);
double xireal(double rval);
int pktoxifree();

