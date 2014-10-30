//the linear xi(r) is needed for the pair-weighting correction to vin(r).

#include "cosmoxi2d.h"
#include "linearpk.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "besselroots.h"

int getxilini(double intkmax, double rmin, double rmax, double deltar, double bias);
int printxilin(char *outfname);
double xilin(double rval);
int pktoxilinfree();
