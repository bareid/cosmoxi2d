#include "xistreamLPT.h"

static double howmanysig = 10.; //this specifies the integration limits -- how many sigma over which to do the integral.
static double sqrt2pi;
static int ns;
static int ellmax; //this is the ellmax used to allocate xis.
static double *s; //vector of redshift space coordinates at which to evaluate xis.
static double **xis;

static xistreamsetup = 0; //switch to 1 after allocating integration workspace

static size_t limitint = 20000;
static double epsabs = 1.0e-7;
static double epsrel = 1.0e-7;
static gsl_integration_workspace *workspaceint;
static gsl_integration_workspace *workspaceintmu;

static gsl_spline **splinexistreamLPT;
static gsl_interp_accel **accxistreamLPT;

int streamLPTsetup() {
  if(xistreamsetup == 1) {
    return 0;
    }

  workspaceint = gsl_integration_workspace_alloc(limitint);
  workspaceintmu = gsl_integration_workspace_alloc(limitint);
  xistreamsetup=1;

  assert(ellmax%2 == 0);
  assert(ns > 0);
  s = (double *) malloc(sizeof(double)*ns);
  xis = malloc2ddouble(ellmax/2+1,ns);
  }

int streamLPTfree() {
  if(xistreamsetup == 1) {
    gsl_integration_workspace_free(workspaceint);
    gsl_integration_workspace_free(workspaceintmu);
    if(xis) {
      malloc2dfree(xis);
      }
    if(s) {
      free(s);
      }
    xistreamsetup = 0;
    }
  return 0;
  }

//los integration stuff.
typedef struct {
double rsig;
double rpi;
double iso1d; //this is added to sig2, independent of mu.
} losint_params;

typedef struct {
double sval;
int ell;
double iso1d; //this is added to sig2, independent of mu.
double alphaperp;
double alphapar;
} muint_params;

static double losint_integrand(double y, void *params) {
  losint_params *p = (losint_params *) params;
  double rsig = p->rsig;
  double rpi = p->rpi;
  double iso1d = p->iso1d;

  double myr = sqrt(rsig*rsig + y*y);
  double mymu = y/myr; //y can be positive or negative, and the sign of the infall velocity los component mymu*vinfall is correct.
  double mymusqr = mymu*mymu;
  if(isnan(mymusqr)) {
    printf("isnan: %e %e %e %e %e\n",y,myr,rsig,mymu,mymusqr);
    }
  double sig2 = iso1d;
  double sqrtsig2 = sqrt(sig2);
  double exparg = (rpi - y)*(rpi - y)/sig2;
  //no need to integrate the constant 1.+\xi, because the dispersion is independent of scale.
  return (xiLPTred(myr,mymu)*exp(-exparg*0.5)/sqrt2pi/sqrtsig2);
  }

static double losint(void *params) {
  double losint_result, losint_err;
  gsl_function AA;
  int status;
  AA.function = &losint_integrand;
  AA.params = params;
  losint_params *p = (losint_params *) params;
  double iso1d = p->iso1d;
  double mysig = sqrt(iso1d);
  double rpi = p->rpi;
  status = gsl_integration_qag(&AA,rpi - howmanysig*mysig, rpi + howmanysig*mysig, epsabs, epsrel, limitint, GSL_INTEG_GAUSS31, workspaceint, &losint_result, &losint_err);
  return losint_result;
  }

static double muint_integrand(double muval, void *params) {
  muint_params *p = (muint_params *) params;
  double sval = p->sval;
  int ell = p->ell;
  double iso1d = p->iso1d;
  double alphaperp = p->alphaperp;
  double alphapar = p->alphapar;

  double rpi = sval*muval*alphapar;
  double rsig = sval*sqrt(1.-muval*muval)*alphaperp;
  losint_params myq = {rsig, rpi, iso1d};
  return((losint((void *) &myq))*legendrep(ell,muval)*(2.*((double) ell)+1.));
  }

static double muint(void *params) {
  double muint_result, muint_err;
  gsl_function BB;
  int status;
  BB.function = &muint_integrand;
  BB.params = params;

  status = gsl_integration_qag(&BB,0.,1., epsabs, epsrel, limitint, GSL_INTEG_GAUSS15, workspaceintmu, &muint_result, &muint_err);
  if(status != 0) {
    //printf("int error: %s\n",gsl_strerror(status));
    }
  return muint_result;
  }

//end los integration stuff.
// fvel is the parameter we want to fit for (fvel = 0.74429 for our a=0.6452 mocks).
int dostreamxiLPT(double smin, double smax, double deltas, double isotropicdisp1d, int ellmaxin, double alphaperp, double alphapar) {
  sqrt2pi = sqrt(2.*M_PI);
  ns = floor((smax - smin)/deltas) + 2;
  if(( smin+deltas*((float) ns) - smax) > 0.99) {
    //ns = ns - 1;
    }
  if(ellmaxin != ellmax) { //deallocate, then reallocate.
    streamfree();
    }
  ellmax = ellmaxin;
  streamLPTsetup();
  gsl_set_error_handler_off();

  int sloop,ello;
  double sval;
  for(sloop=0;sloop<ns;sloop++) {
    sval = smin+((float) (sloop))*deltas;
    s[sloop] = sval;
    for(ello=0;ello<=ellmax;ello+=2) {
      muint_params mypmu = {sval,ello,isotropicdisp1d,alphaperp,alphapar};
      xis[ello/2][sloop] = muint((void *) &mypmu);
      }
    }

  //initialize spline.
  splinexistreamLPT = (gsl_spline **) malloc((ellmax/2+1)*sizeof(gsl_spline *));
  accxistreamLPT = (gsl_interp_accel **) malloc((ellmax/2+1)*sizeof(gsl_interp_accel *));
  int j;
  for(j=0;j<=ellmax/2;j++) {
    splinexistreamLPT[j] = gsl_spline_alloc(gsl_interp_linear,ns);
    accxistreamLPT[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinexistreamLPT[j], s, xis[j], ns);
    }
  return 0;
  }

double xistreamLPTeval(double sval, int ell) {
  //no tests for if we're outside spline boundaries!!
  if(sval < s[0] || sval > s[ns-1]) {
    printf("xistreamLPTeval outside computed bounds %e %e %e\n",sval,s[0],s[ns-1]);
    exit(1);
    }
  return(gsl_spline_eval(splinexistreamLPT[ell/2],sval,accxistreamLPT[ell/2]));
  }

int printxisLPT(char *outfname, double bias, double fvel, double iso1d, double alphaperp, double alphapar) {
  int sloop, ello;
  FILE *ofp;
  ofp = open_file_write(outfname);
  fprintf(ofp,"b = %f\n",bias);
  fprintf(ofp,"f = %f\n",fvel);
  fprintf(ofp,"iso1d = %f\n",iso1d);
  fprintf(ofp,"alphaperp = %f\n",alphaperp);
  fprintf(ofp,"alphapar = %f\n",alphapar);
  for(sloop=0;sloop<ns;sloop++) {
    fprintf(ofp,"%e ",s[sloop]);
    for(ello=0;ello<=ellmax;ello+=2) {
      fprintf(ofp,"%e ",xis[ello/2][sloop]);
      } //end ello
    fprintf(ofp,"\n");
    } //end sloop
  fclose(ofp);
  return 0;
  } //end printxis
