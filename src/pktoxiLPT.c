#include "pktoxiLPT.h"

static int xisetup = 0; //switch to 1 after allocating integration workspace and filling besselroots.

static size_t limitkint=20000;
static double epsabskint=1.0e-9; //might need to adjust this later.
static double epsrelkint=1.e-6;
static gsl_integration_workspace *workspacekint;

static double *r;  //these are really in redshift space, but for copy/paste, we'll just leave them labeled by r instead ofs.
static double **xi; //xi contains all the ells.
static int nr;
static int xiLPTellmax;

static gsl_interp_accel **accxi;
static gsl_spline **splinexi;


typedef struct {
double rval;
int ell;
} kint_params;

int doxiLPTsetup(int ellmax) {
  if(xisetup == 1) {
    return 0;
    }

  fillbesselrootsxiLPT(ellmax); //done in besselroots.h
  workspacekint =  gsl_integration_workspace_alloc (limitkint);
  r = malloc(sizeof(double)*nr);
  xi = malloc2ddouble(ellmax/2+1,nr);

  //splines allocated below.
  xisetup = 1;
  }

int pktoxiLPTfree() {
  int j;
  if(workspacekint) {
    gsl_integration_workspace_free(workspacekint);
    }
  if(r) {
    free(r);
    }
  if(xi) {
    malloc2dfree(xi);
    }
  if(splinexi[0]) {
    for(j=0;j<=xiLPTellmax/2;j++) {
      gsl_spline_free(splinexi[j]);
      gsl_interp_accel_free(accxi[j]);
      }
    free(splinexi);
    free(accxi); 
    }
  }

//labelled static so only available
static double kint_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  int ell = paramsk->ell;
  return (kval*kval*pkredLPT(kval,ell)*gsl_sf_bessel_jl(ell,kval*rval));
  }

//input the max k to integrate to, independent of rval here.

static double kint(double intkmax, void *params) {
  kint_params *paramsk = (kint_params *) params;
  double rval = paramsk->rval;
  int ell = paramsk->ell;

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
  if(getbesselrootxiLPT(whichroot,ell) < 1.0e-4) {
    whichroot += 1;
    }
  root = getbesselrootxiLPT(whichroot,ell)/rval;
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
    root = getbesselrootxiLPT(whichroot,ell)/rval;
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

  double ellprefac;
  assert(ell%2==0);
  if(ell%4==0) {
    ellprefac = 1.;
    }
  else {
    ellprefac = -1.;
    }

  return (myinttot*ellprefac/2./M_PI/M_PI);
  }

//intkmax is max k for which pkLPT(k) > 0.
int getxiLPTmoments(double intkmax, double smin, double smax, double deltas, int ellmax) {
  xiLPTellmax = ellmax;
  assert(xiLPTellmax <= 6);  //if i want to go higher, I have to define higher legendre polys in legendrep(ell,x) in misc.c
  nr = floor((smax - smin)/deltas) + 2;
  if(( smin+deltas*((float) nr) - smax) > 0.99) {
    nr = nr - 1;
    }
  doxiLPTsetup(ellmax);
  gsl_set_error_handler_off();
  int rloop,j;
  double rval;
  for(rloop=0;rloop<nr;rloop++) {
    rval = smin + ((double) rloop)*deltas;
    r[rloop] = rval;
    for(j=0;j<=ellmax/2;j++) {
      kint_params myp = {rval, j*2};
      xi[j][rloop] = kint(intkmax,(void *) &myp);
      }
    }
  //setup spline.

  splinexi = (gsl_spline **) malloc((ellmax/2+1)*sizeof(gsl_spline *));
  accxi = (gsl_interp_accel **) malloc((ellmax/2+1)*sizeof(gsl_interp_accel *));

  for(j=0;j<=ellmax/2;j++) {
    splinexi[j] = gsl_spline_alloc(gsl_interp_linear,nr);
    accxi[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinexi[j], r, xi[j], nr);    
    }
  return 0;
  }

double xiLPTredmoments(double sval, int ell) {
  //hopefully this one never goes outside.
  assert(sval >= r[0] && sval <= r[nr-1]);
  return (gsl_spline_eval(splinexi[ell/2],sval,accxi[ell/2]));
  }

int printxiLPTmoments(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i,j;
  for(i=0;i<nr;i++) {
    fprintf(ofp,"%e ",r[i]);
    for(j=0;j<=xiLPTellmax/2;j++) {
      fprintf(ofp,"%e ",xi[j][i]);
      }
    fprintf(ofp,"\n");
    }
  fclose(ofp);
  return 0;
  }

double xiLPTred(double sval, double mu) {
  double xival = 0.;
  int j,ell;
  for(j=0;j<=xiLPTellmax/2;j++) {
    ell = 2*j;
    xival += xiLPTredmoments(sval,ell)*legendrep(ell,mu);
    }
  return xival;
  }
