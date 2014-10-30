#include "pktoxi.h"

static int xisetup = 0; //switch to 1 after allocating integration workspace and filling besselroots.

static double sqrtpi;

static size_t limitkint=20000;
static double epsabskint=1.0e-9; //might need to adjust this later.
static double epsrelkint=1.e-6;
static gsl_integration_workspace *workspacekint;

static double *r;
static double *xi;
static int nr;

static gsl_interp_accel *accxi;
static gsl_spline *splinexi;


typedef struct {
double rval;
} kint_params;

int doxisetup() {
  if(xisetup == 1) {
    return 0;
    }

  fillbesselroots(); //done in besselroots.h
  sqrtpi = sqrt(M_PI);
  workspacekint =  gsl_integration_workspace_alloc (limitkint);
  r = malloc(sizeof(double)*nr);
  xi = malloc(sizeof(double)*nr);
  accxi = gsl_interp_accel_alloc ();
  splinexi = gsl_spline_alloc (gsl_interp_linear, nr);
  xisetup = 1;
  }

int pktoxifree() {
  if(workspacekint) {
    gsl_integration_workspace_free(workspacekint);
    }
  if(r) {
    free(r);
    }
  if(xi) {
    free(xi);
    }
  }

//labelled static so only available
static double kint_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  return (kval*kval*pkLPT(kval)*gsl_sf_bessel_j0(kval*rval));
  }

//input the max k to integrate to, independent of rval here.

static double kint(double intkmax, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;

  double kInt_err;
  gsl_function KK;
  int kstatus;
  KK.function = &kint_integrand;
  KK.params = params;

  double kmin, kmax;
  double myint,myinttot,myintpos,myintneg;
  double myintposerr, myintnegerr;
  double root;
  int whichroot=0;

  kmin = 0.0;
  if(getbesselroot(whichroot,0) < 1.0e-4) {
    whichroot += 1;
    }
  root = getbesselroot(whichroot,0)/rval;
  kmax = root;
  whichroot += 1;

  myinttot = 0.;
  myintpos = 0.;
  myintneg = 0.;
  myintposerr = 0.;
  myintnegerr = 0.;

    while(kmin < intkmax) { //get all the way to kmax.
    kstatus = gsl_integration_qag(&KK,kmin,kmax,epsabskint,epsrelkint,limitkint,GSL_INTEG_GAUSS31,workspacekint, &myint, &kInt_err);
    myinttot += myint;
    if(myint > 0) {
      myintpos += myint;
      myintposerr += kInt_err;
      }
    else {
      myintneg += myint;
      myintnegerr += kInt_err;
      }
    kmin = kmax;
    root = getbesselroot(whichroot,0)/rval;
    kmax = root;
    whichroot += 1;
    assert(kmin < kmax);
    }  //end while

  if(fabs(myinttot - myintpos - myintneg) > epsabskint && fabs(myinttot - myintpos - myintneg) > fabs(epsrelkint*myinttot) ) {
    printf("CONVERGENCE PROB??? %e %e %e %e %e\n",myinttot,myintpos,myintneg,epsabskint,epsrelkint);
    exit(1);
    }

  if((myinttot < 0. && fabs(myinttot) < 1.0e-5*fabs(myintneg)) || (myinttot > 0. && myinttot < 1.0e-5*myintpos)) {
    printf("CONVERGENCE PROB2??? %e %e %e %e %e\n",myinttot,myintpos,myintneg,1.0e-5,1.0e-5);
    }

  return (myinttot/2./M_PI/M_PI);
  }

//intkmax is max k for which pkLPT(k) > 0.
int getxireal(double intkmax, double rmin, double rmax, double deltar) {
  nr = floor((rmax - rmin)/deltar) + 2;
  if(( rmin+deltar*((float) nr) - rmax) > 0.99) {
    nr = nr - 1;
    }
  doxisetup();
  gsl_set_error_handler_off();
  int rloop;
  double rval;
  for(rloop=0;rloop<nr;rloop++) {
    rval = rmin + ((double) rloop)*deltar;
    r[rloop] = rval;
    kint_params myp = {rval};
    xi[rloop] = kint(intkmax,(void *) &myp);
    }
  //setup spline.
  gsl_spline_init(splinexi,r,xi,nr);
  return 0;
  }

int printxireal(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i;
  for(i=0;i<nr;i++) {
    fprintf(ofp,"%e %e\n",r[i],xi[i]);
    }
  fclose(ofp);
  return 0;
  }

double xireal(double rval) {
  if(rval < r[0]) {
    return xi[0];
    }
  if(rval > r[nr-1]) {
    assert(rval < 1.05*r[nr-1]);
    double slope = (xi[nr-1] - xi[nr - 2])/(r[nr-1] - r[nr-2]);
    return(xi[nr-1] + slope*(rval - r[nr-1]));
    }
  return (gsl_spline_eval(splinexi,rval,accxi));
  }
