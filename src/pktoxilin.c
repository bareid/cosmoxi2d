#include "pktoxilin.h"

static int xilinsetup = 0; //switch to 1 after allocating integration workspace and filling besselroots.

static double sqrtpi;

static size_t limitkint=20000;
static double epsabskint=1.0e-9; //might need to adjust this later.
static double epsrelkint=1.e-6;
static gsl_integration_workspace *workspacekint;

static double *rlin;
static double *xilini;
static int nrlin;

static gsl_interp_accel *accxilin;
static gsl_spline *splinexilin;

static const double sigmasqr = 1.0e-3;  //helps the integral converge at high k.

typedef struct {
double rval;
} kint_params;

int doxilinsetup() {
  if(xilinsetup == 1) {
    return 0;
    }

  fillbesselroots(); //done in besselroots.h
  sqrtpi = sqrt(M_PI);
  workspacekint =  gsl_integration_workspace_alloc (limitkint);
  rlin = malloc(sizeof(double)*nrlin);
  xilini = malloc(sizeof(double)*nrlin);
  accxilin = gsl_interp_accel_alloc ();
  splinexilin = gsl_spline_alloc (gsl_interp_linear, nrlin);
  xilinsetup = 1;
  }

int pktoxilinfree() {
  if(workspacekint) {
    gsl_integration_workspace_free(workspacekint);
    }
  if(rlin) {
    free(rlin);
    }
  if(xilini) {
    free(xilini);
    }
  }

//labelled static so only available
static double kint_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  return (kval*kval*pk(kval)*gsl_sf_bessel_j0(kval*rval)*exp(-kval*kval*sigmasqr)); //this is the linear power spectrum.
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

  while(kmax < intkmax) {
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

int getxilini(double intkmax, double rmin, double rmax, double deltar, double bias) {
  assert(intkmax > 5.);  //this was 10. in file I copied from.
  nrlin = floor((rmax - rmin)/deltar) + 2;
  if(( rmin+deltar*((float) nrlin) - rmax) > 0.99) {
    nrlin = nrlin - 1;
    }
  doxilinsetup();
  gsl_set_error_handler_off();
  int rloop;
  double rval;
  for(rloop=0;rloop<nrlin;rloop++) {
    rval = rmin + ((double) rloop)*deltar;
    rlin[rloop] = rval;
    kint_params myp = {rval};
    xilini[rloop] = kint(intkmax, (void *) &myp)*bias*bias;
    }
  //setup spline.
  gsl_spline_init(splinexilin,rlin,xilini,nrlin);
  return 0;
  }

int printxilin(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i;
  for(i=0;i<nrlin;i++) {
    fprintf(ofp,"%e %e\n",rlin[i],xilini[i]);
    }
  fclose(ofp);
  return 0;
  }

double xilin(double rval) {
  if(rval < rlin[0]) {
    return xilini[0];
    }
  if(rval > rlin[nrlin-1]) {
    assert(rval < 1.05*rlin[nrlin-1]);
    double slope = (xilini[nrlin-1] - xilini[nrlin - 2])/(rlin[nrlin-1] - rlin[nrlin-2]);
    return(xilini[nrlin-1] + slope*(rval - rlin[nrlin-1]));
    }
/*
  if(rval > rlin[nrlin-1]) {
    printf("warning r > r[nr-1], returning xi(%e) = %e, xi(r[0]=%e) = %e\n",rval,gsl_spline_eval(splinexilin,rval,accxilin),rlin[nrlin-1],gsl_spline_eval(splinexilin,rlin[nrlin-1],accxilin));
    }
*/
  return (gsl_spline_eval(splinexilin,rval,accxilin));
  }
