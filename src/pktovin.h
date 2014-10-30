#include "cosmoxi2d.h"
#include "linearpk.h"
#include "ptraw.h"
#include "pktoxilin.h" //may want to replace with xiLPT later.
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "besselroots.h"

//this now computed vinlin by default as well.  return to this later and see if we can replace with the nonlinear version.
int getvinreal(double intkmax, double rmin, double rmax, double deltar);
int getvinlinreal(double intkmax, double rmin, double rmax, double deltar, double b1L);
int printvinreal(char *outfname);
double vinreal(double rval);
int printvinlinreal(char *outfname);
double vinlinreal(double rval);
int pktovinfree();
