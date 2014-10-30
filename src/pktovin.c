#include "pktovin.h"

static int vinsetup = 0;

static double sqrtpi;

static size_t limitkint=20000;
static double epsabskint=1.0e-9; //might need to adjust this later.
static double epsrelkint=1.e-6;
static gsl_integration_workspace *workspacekint;

static double *r;
static double *rlin;
static double *vin;
static double *vinlin;
static int nr=0;

static gsl_interp_accel *accvin;
static gsl_spline *splinevin;

static gsl_interp_accel *accvinlin;
static gsl_spline *splinevinlin;


typedef struct {
double rval;
int nlorlin; //nl = 0, lin = 1.
} kint_params;

int dovinsetup() {
  if(vinsetup == 1) {
    return 0;
    }

  fillbesselroots();
  sqrtpi = sqrt(M_PI);
  workspacekint =  gsl_integration_workspace_alloc (limitkint);
  r = malloc(sizeof(double)*nr);
  rlin = malloc(sizeof(double)*nr);
  vin = malloc(sizeof(double)*nr);
  vinlin = malloc(sizeof(double)*nr);
  accvin = gsl_interp_accel_alloc ();
  splinevin = gsl_spline_alloc (gsl_interp_linear, nr);
  accvinlin = gsl_interp_accel_alloc ();
  splinevinlin = gsl_spline_alloc (gsl_interp_linear, nr);
  vinsetup = 1;
  }

int pktovinfree() {
  if(workspacekint) {
    gsl_integration_workspace_free(workspacekint);
    }
  if(r) {
    free(r);
    }
  if(vin) {
    free(vin);
    }
  if(vinlin) {
    free(vinlin);
    }
  }

static double kint_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  int nlorlin = paramsk->nlorlin;
  if(nlorlin == 0) {
    return(vinSPTkint(kval)*gsl_sf_bessel_jl(1,kval*rval));
    }
  if(nlorlin == 1) {
    return(pk(kval)*kval*gsl_sf_bessel_jl(1,kval*rval));
    }
  //only nlorlin = 0,1 defined.
  exit(1);
  }

//input the max k to integrate to, independent of rval here.

static double kint(double intkmax, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  int nlorlin = paramsk->nlorlin;

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
    assert(kmin < kmax); //if you get to end of bessel roots, kmax gets assigned 0.
    }  //end while

  if(fabs(myinttot - myintpos - myintneg) > epsabskint && fabs(myinttot - myintpos - myintneg) > fabs(epsrelkint*myinttot) ) {
    printf("CONVERGENCE PROB??? %e %e %e %e %e\n",myinttot,myintpos,myintneg,epsabskint,epsrelkint);
    exit(1);
    }

  if((myinttot < 0. && fabs(myinttot) < 1.0e-5*fabs(myintneg)) || (myinttot > 0. && myinttot < 1.0e-5*myintpos)) {
    printf("CONVERGENCE PROB2??? %e %e %e %e %e\n",myinttot,myintpos,myintneg,1.0e-5,1.0e-5);
    }

  return (myinttot/M_PI/M_PI);
  }

//intkmax is max k for which pkLPT(k) > 0.
int getvinreal(double intkmax, double rmin, double rmax, double deltar) {
  int nrtest = floor((rmax - rmin)/deltar) + 2;
  if(( rmin+deltar*((float) nrtest) - rmax) > 0.99) {
    nrtest = nrtest - 1;
    }
  if(nr == 0) {
    nr = nrtest;
    }
  else {
    assert(nr == nrtest);
    }
  dovinsetup();
  gsl_set_error_handler_off();
  int rloop;
  double rval;
  for(rloop=0;rloop<nr;rloop++) {
    rval = rmin + ((double) rloop)*deltar;
    r[rloop] = rval;
    kint_params myp = {rval,0};
    vin[rloop] = kint(intkmax,(void *) &myp)/(1.+xilin(rval));
    }
  //setup spline.
  gsl_spline_init(splinevin,r,vin,nr);
  return 0;
  }

int getvinlinreal(double intkmax, double rmin, double rmax, double deltar, double b1L) {
  int nrtest = floor((rmax - rmin)/deltar) + 2;
  if(( rmin+deltar*((float) nrtest) - rmax) > 0.99) {
    nrtest = nrtest - 1;
    }
  if(nr == 0) {
    nr = nrtest;
    }
  else {
    assert(nr == nrtest);
    }

  dovinsetup();
  gsl_set_error_handler_off();
  int rloop;
  double rval;
  for(rloop=0;rloop<nr;rloop++) {
    rval = rmin + ((double) rloop)*deltar;
    rlin[rloop] = rval;
    kint_params myplin = {rval,1};
    vinlin[rloop] = (1.+b1L)*kint(intkmax,(void *) &myplin);
    }
  //setup spline.
  gsl_spline_init(splinevinlin,rlin,vinlin,nr);
  return 0;
  }

int printvinreal(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i;
  for(i=0;i<nr;i++) {
    fprintf(ofp,"%e %e\n",r[i],vin[i]);
    }
  fclose(ofp);
  return 0;
  }

int printvinlinreal(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i;
  for(i=0;i<nr;i++) {
    fprintf(ofp,"%e %e\n",r[i],vinlin[i]);
    }
  fclose(ofp);
  return 0;
  }

double vinreal(double rval) {
  if(rval < r[0]) {
    return vin[0];
    }
  if(rval > r[nr-1]) {
    double slope = (vin[nr-1] - vin[nr - 2])/(r[nr-1] - r[nr-2]);
    assert(rval < r[nr-1]*1.05);
    return(vin[nr-1] + slope*(rval - r[nr-1]));
    }
  return (gsl_spline_eval(splinevin,rval,accvin));
  }

double vinlinreal(double rval) {
  if(rval < rlin[0]) {
    return vin[0];
    }
  if(rval > rlin[nr-1]) {
    double slope = (vinlin[nr-1] - vinlin[nr - 2])/(rlin[nr-1] - rlin[nr-2]);
    assert(rval < rlin[nr-1]*1.05);
    return(vinlin[nr-1] + slope*(rval - rlin[nr-1]));
    }
  return (gsl_spline_eval(splinevinlin,rval,accvinlin));
  }
