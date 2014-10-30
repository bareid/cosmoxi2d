#include "cosmoxi2d.h"
#include "linearpk.h"
#include "ptraw.h"
//sigsqr does not need lienar xi, cancels out for next-to-leading order.  does need linear theory velocity though.
//#include "pktoxilin.h" //may want to replace with xiLPT later.
#include "pktovin.h"

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "besselroots.h"

int pktosigsqrfree();
int getsigsqrreal(double intkmax, double rmin, double rmax, double deltar, double b1L);
int printsigsqr(char *outfname);
double sigsqrperpreal(double rval);
double sigsqrparreal(double rval);
double sigsqrtot(double musqr, double rval);
