//from SPTveldisp/iset1.c and iset2.c
#include "pktosig.h"

static int sigsqrsetup = 0;

static double sqrtpi;

static size_t limitkint=20000;
static double epsabskint=1.0e-9; //might need to adjust this later.
static double epsrelkint=1.e-6;
static gsl_integration_workspace *workspacekint;

static double *r;
static double *sigsqrperp;
static double *sigsqrpar;
static int nr;

static gsl_interp_accel *accsigsqrperp;
static gsl_spline *splinesigsqrperp;

static gsl_interp_accel *accsigsqrpar;
static gsl_spline *splinesigsqrpar;

typedef struct {
int whichbessel;
int whichkint; //which integrand
double rval;
} kint_params;


int dosigsqrsetup() {
  if(sigsqrsetup == 1) {
    return 0;
    }

  fillbesselroots();
  fillbesselrootsA13();
  sqrtpi = sqrt(M_PI);
  workspacekint =  gsl_integration_workspace_alloc (limitkint);
  r = malloc(sizeof(double)*nr);
  sigsqrperp = malloc(sizeof(double)*nr);
  sigsqrpar = malloc(sizeof(double)*nr);
  accsigsqrperp = gsl_interp_accel_alloc ();
  accsigsqrpar = gsl_interp_accel_alloc ();
  splinesigsqrperp = gsl_spline_alloc (gsl_interp_linear, nr);
  splinesigsqrpar = gsl_spline_alloc (gsl_interp_linear, nr);
  sigsqrsetup = 1;
  }

int pktosigsqrfree() {
  if(workspacekint) {
    gsl_integration_workspace_free(workspacekint);
    }
  if(r) {
    free(r);
    }
  if(sigsqrperp) {
    free(sigsqrperp);
    }
  if(sigsqrpar) {
    free(sigsqrpar);
    }
  }

double kint_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;

  int whichbessel = paramsk->whichbessel;
  int whichkint = paramsk->whichkint;
  double rval = paramsk->rval; 
  return(mybessel(whichbessel,kval*rval)*sigsqrSPTkint(kval,whichkint));
  }


double kintnowig(double intkmax) {


  int whichbessel = -1;
  int whichkint = 0;
  kint_params mypx = {whichbessel,whichkint,0.0};

  double kInt_err,myint;
  gsl_function KK;
  int kstatus;
  KK.function = &kint_integrand;
  KK.params = &mypx;

  kstatus = gsl_integration_qag(&KK,0.0,intkmax,epsabskint,epsrelkint,limitkint,GSL_INTEG_GAUSS31,workspacekint, &myint, &kInt_err);

  //does not need to be super precise because we're going to marginalize over an additional sigma2.
  return(myint/M_PI/M_PI/3.);
  }

double kint(double intkmax, void *params) {
  kint_params *paramsk = (kint_params *) params;

  int whichbessel = paramsk->whichbessel;
  int whichkint = paramsk->whichkint;
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

  int ell;
  if(whichbessel == 0) {
    ell=0;
    }
  else {
    if(whichbessel == 1) {
      ell=1;
      }
    else {
      exit(1);
      }
    }

  kmin = 0.0;
  root = getbesselroot(whichroot,ell)/rval;
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
    root = getbesselroot(whichroot,ell)/rval;
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
  return(myinttot/M_PI/M_PI);
  }

double kintA13_integrand(double kval, void *params) {
  kint_params *paramsk = (kint_params *) params;
  int whichbessel = paramsk->whichbessel;
  int kopt = paramsk->whichkint;
  double rval = paramsk->rval;
  switch(kopt) {
    case(0):
  return(pk(kval)*mybesselA13(whichbessel,kval*rval));
    case(1):
  return(kval*pk(kval)*mybesselA13(whichbessel,kval*rval));
    case(2):
  return(gsl_pow_2(kval)*pk(kval)*mybesselA13(whichbessel,kval*rval));
    default:
  exit(1);
    }
  }

double kintA13(double intkmax, void *params) {
  kint_params *paramsk = (kint_params *) params;
  int whichbessel = paramsk->whichbessel;
  int kopt = paramsk->whichkint;
  double rval = paramsk->rval;

  double kInt_err;
  gsl_function KK;
  int kstatus;
  KK.function = &kintA13_integrand;
  KK.params = params;

  double kmin, kmax;
  double myint,myinttot,myintpos,myintneg;
  double myintposerr, myintnegerr;
  double root;
  int whichroot=0;

  kmin = 0.0;
  root = getbesselrootA13(whichroot,whichbessel)/rval;
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
    root = getbesselrootA13(whichroot,whichbessel)/rval;
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

  return (myinttot);
  }

//since b is already put into sigsqr, need to put it in here as well.
int getA13(double intkmax, double rval, double *a13perp, double *a13par) {
  int koptlist[7];  //k to the power 0, 1, 2 in the integrand.
  int whichbessellist[7]; //which oscillatory function.

  double xiresult[7];
  //0: k P(k) j1(kr)
  //these four are needed for line 3 of eqn 30; can be reexpressed in terms of \xi(r), linear velocity correlation functions, and v12(r)/r, but i'll just brute force here.  see notes sigma_12^2(r) pg 51
  //1: k^2 P(k) j0(kr)
  //2: k^2 P(k) j1(kr)/kr
  //3: k^0 P(k) j0(kr)
  //4: k^0 P(k) j1(kr)/kr
  //5: k P(k) I13(r)
  //6: k P(k) I31(r)

  koptlist[0] = 1;
  koptlist[5] = 1;
  koptlist[6] = 1;

  koptlist[1] = 2;
  koptlist[2] = 2;

  koptlist[3] = 0;
  koptlist[4] = 0;

  whichbessellist[0] = 2;
  whichbessellist[1] = 0;
  whichbessellist[2] = 1;
  whichbessellist[3] = 0;
  whichbessellist[4] = 1;
  whichbessellist[5] = 3;
  whichbessellist[6] = 4;

  int bigint;
  for(bigint=0;bigint<7;bigint++) {
    kint_params myp2 = {whichbessellist[bigint],koptlist[bigint],rval};
    xiresult[bigint] = kintA13(intkmax, (void *) &myp2);
    }

  (*a13perp) = -(0.5*xiresult[2]*xiresult[4]-2./7.*xiresult[6]*xiresult[6])/gsl_pow_4(M_PI);
  (*a13par) = -(-5./14.*xiresult[0]*xiresult[0]+0.5*(xiresult[1]-2.*xiresult[2])*(xiresult[3]-2.*xiresult[4])-1./7.*xiresult[5]*xiresult[5]-2./7.*xiresult[6]*xiresult[6])/gsl_pow_4(M_PI);
  return 0;
  }

int getsigsqrreal(double intkmax, double rmin, double rmax, double deltar, double b1L) {
  nr = floor((rmax - rmin)/deltar) + 2;
  if(( rmin+deltar*((float) nr) - rmax) > 0.99) {
    nr = nr - 1;
    }
  dosigsqrsetup();
  gsl_set_error_handler_off();
  int rloop;
  double rval;
  double myi[3][2]; //first index is which k integral, second is which bessel
  double a13perp,a13par;
  double sig2infty = kintnowig(intkmax);
  #ifdef VERBOSE
  printf("disp at infty = %e\n",sig2infty);
  #endif
  int whichkint, whichbessel;

  for(rloop=0;rloop<nr;rloop++) {
    rval = rmin + ((double) rloop)*deltar;
    r[rloop] = rval;
    getA13(intkmax, rval, &a13perp, &a13par);
    for(whichbessel=0;whichbessel<=1;whichbessel++) {
      for(whichkint=0;whichkint<=2;whichkint++) {
        if(whichkint == 2 && whichbessel == 1) {
          myi[whichkint][whichbessel] = 0.;
          continue;
          }
        kint_params myp1 = {whichbessel,whichkint,rval};
        myi[whichkint][whichbessel] = kint(intkmax, &myp1);
        }
      }
    sigsqrperp[rloop] = (sig2infty - myi[0][1]) + myi[2][0] + myi[1][1] + 2.*(1.+b1L)*a13perp;
    sigsqrpar[rloop] = (sig2infty - (myi[0][0] - 2.*myi[0][1])) + myi[2][0] + (myi[1][0] - 2.*myi[1][1]) + 2.*(1.+b1L)*a13par - 0.5*gsl_pow_2(vinlinreal(rval));
    //these can be negative in smallest rbins.
    //Lado solution is to set them to a small positive value:
    if(sigsqrperp[rloop] <= 0.) {
      sigsqrperp[rloop] = 0.1;
      }
    if(sigsqrpar[rloop] <= 0.) {
      sigsqrpar[rloop] = 0.1;
      }
//    printf("%e %e %e\n",rval,sigsqrperp[rloop],sigsqrpar[rloop]);
    }
  gsl_spline_init(splinesigsqrpar,r,sigsqrpar,nr);
  gsl_spline_init(splinesigsqrperp,r,sigsqrperp,nr);
  return 0;
  }

int printsigsqr(char *outfname) {
  FILE *ofp = open_file_write(outfname);
  int i;
  for(i=0;i<nr;i++) {
    fprintf(ofp,"%e %e %e\n",r[i],sigsqrperp[i],sigsqrpar[i]);
    }
  fclose(ofp);
  return 0;
  }

double sigsqrperpreal(double rval) {
  if(rval < r[0]) {
    return sigsqrperp[0];
    }
  if(rval > r[nr-1]) { 
    double slope = (sigsqrperp[nr-1] - sigsqrperp[nr - 2])/(r[nr-1] - r[nr-2]);
    assert(rval < r[nr-1]*1.05);
    return(sigsqrperp[nr-1] + slope*(rval - r[nr-1]));
    }
  //else, within the allowed splining region.
  return (gsl_spline_eval(splinesigsqrperp,rval,accsigsqrperp));
  }

double sigsqrparreal(double rval) {
  if(rval < r[0]) {
    return sigsqrpar[0];
    }
  if(rval > r[nr-1]) { 
    double slope = (sigsqrpar[nr-1] - sigsqrpar[nr - 2])/(r[nr-1] - r[nr-2]);
    assert(rval < r[nr-1]*1.05);
    return(sigsqrpar[nr-1] + slope*(rval - r[nr-1]));
    }
  //else, within the allowed splining region.
  return (gsl_spline_eval(splinesigsqrpar,rval,accsigsqrpar));
  }

//PT predicts negative velocity dispersion on really small scales -- we need to suppress that to avoid nan's.  Just set to 0.
double sigsqrtot(double musqr, double rval) {
  if(!(musqr >= 0. && musqr <= 1.)) {
    printf("weird musqr %e\n",musqr);
    }
  assert(musqr >= 0. && musqr <= 1.);
  double mysigsqrtot = (1.-musqr)*sigsqrperpreal(rval) + musqr*sigsqrparreal(rval);
  if(mysigsqrtot < 0.) {
    assert(0==1);  //should have solved this problem already!
    return(0.);
    }
  return ((1.-musqr)*sigsqrperpreal(rval) + musqr*sigsqrparreal(rval));
  }
