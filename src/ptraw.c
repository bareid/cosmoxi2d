#include "ptraw.h"

static int ptsetup = 0; //switch to 1 after allocating integration workspaces.

static int lenkqr;
static double *klistqr;
static double **qlist;  //2d array for the Q integral results.
static double **rlist;  //2d array for the R integral results.

//splines for q,r
static gsl_spline **splineqlist;
static gsl_interp_accel **accqlist;
static gsl_spline **splinerlist;
static gsl_interp_accel **accrlist;

static int lenk = 0;  //must be set when generating pkLPTreal or pkLPTred.
static double *klist;  //needs to get assigned somewhere!
static double *pkLPTreal;


static int lenkvel = 0; //this is going to be the fine binning of k, coarser to evaluate PT integrals.
static double *klistvel;
static double *logklistvel;
static double *logvISPT;  //the function of k integrated against j1(kr) to get vin(r).  We're going to interpolate in log
static double *logsigsqrSPT0; //Ptt(k) in perturbation theory.  to be integrated against J(mu^2,r) and 1/3.  
static double *logsigsqrSPT1;  //to be integrated against J(mu^2,r)
static double *logsigsqrSPT2;  //to be integrated against j0(kr) only.
static double b1Lhold = -100; //this is necessary for interpolating vISPT below kvel[0]
static double fvelhold = -100; //this is necessary for interpolating vISPT below kvel[0]
static double nasympvISPT;
static double nasympsigsqrSPT0;
static double nasympsigsqrSPT1;
static double nasympsigsqrSPT2;
//static double *pkLPTred;

static int lenkvelPT = 0; //this is going to be the fine binning of k, coarser to evaluate PT integrals.
static double *klistvelPT;
static double **vinlist; //length 2: P_dt(k) [b  term], Eqn A5 + Eqn A6 [b^2 term] from Reid and White 2011, MNRAS 417, 1913
static double **sig2list; //length 5: P_tt(k), 2 integrals in A12, A14, 15.

static gsl_spline **splinevinlist;
static gsl_interp_accel **accvinlist;
static gsl_spline **splinesig2list;
static gsl_interp_accel **accsig2list;

//integrals in A13 don't need to be stored anywhere in this file holding onto functions of k.

//adding P^LPT(k) in redshift space because it's more accurate at large separations.  Number of elements is ellmax.
static int pkLPTredellmax; //assign this, and use for checks to make sure we don't access outside this LMAX.
static double **pkLPTred; //use same klist as pkLPTreal.
static gsl_spline **splinepkLPTred;
static gsl_interp_accel **accpkLPTred;

//power spectrum spline.
static gsl_interp_accel *accpkLPTreal;
static gsl_spline *splinepkLPTreal;

static gsl_interp_accel *accvISPT;
static gsl_spline *splinevISPT;

static gsl_interp_accel *accsigsqrSPT2;
static gsl_spline *splinesigsqrSPT2;

static gsl_interp_accel *accsigsqrSPT1;
static gsl_spline *splinesigsqrSPT1;

static gsl_interp_accel *accsigsqrSPT0;
static gsl_spline *splinesigsqrSPT0;

/***********************
  NOTES

  1) may want to play with spline type (currently set to linear in k-P(k) space, NOT log-log space!).
  right now we're so densely sampled in k that it probably doesn't matter, but we'll want to tweak it later.

***********************/

//doing one or two-dimensional integrals in this file over P(k*x) P(k*y)
//perform I1 first, then I2.

/********

Begin list of integrals.
1. Q(1,2,...,13), defined on A39.  I1 will be over variable r, and I2 will be over variable x.
********/

size_t limitI1=20000;
double epsabsI1=1.0e-8;
double epsrelI1=1.0e-5;  //adjust later.
static gsl_integration_workspace *workspaceI1;

size_t limitI2=20000;
double epsabsI2=1.0e-8;
double epsrelI2=1.0e-5;  //adjust later.
static gsl_integration_workspace *workspaceI2;

typedef struct {
  double kval;
  double i2val; //value of outer integrand variable.
  int n;
  } i1int_params;

typedef struct {
  double kval;
  int n;
  } i2int_params;

typedef struct {
  double *a; //array with all the coefficients we need.
  double kval;
  double b1L;
  double fvel;
  int ellval;
  } muint_params;

int ptintegralsetup() {
  if(ptsetup == 1) {
    return 0;
    }

  workspaceI1 = gsl_integration_workspace_alloc(limitI1);
  workspaceI2 = gsl_integration_workspace_alloc(limitI2);

  ptsetup = 1;
  return 0;
  }

int ptintegralfree() {
  if(ptsetup == 1) {
    gsl_integration_workspace_free(workspaceI1);
    gsl_integration_workspace_free(workspaceI2);
    ptsetup = 0;
    }
  return 0;
  }

double muint_integrand(double muval, void *params) {
  muint_params *paramsmu = (muint_params *) params;
  double *a = malloc(sizeof(double)*10);
  int j;
  for(j=0;j<10;j++) {
    a[j] = paramsmu->a[j];
    }
  double k = paramsmu->kval;
  double b1L = paramsmu->b1L;
  double fvel = paramsmu->fvel;
  int ell = paramsmu->ellval;
  double musqr = muval*muval;
  return(exp(-(1.+fvel*(fvel+2.)*musqr)*gsl_pow_2(k)*a[0])*(gsl_pow_2(1.+b1L + fvel*musqr)*pk(k) + a[1] + musqr*fvel*(a[2] + fvel*a[3]) + gsl_pow_2(musqr*fvel)*(a[4]+fvel*a[5]+fvel*fvel*a[6]) + gsl_pow_3(musqr*fvel)*(a[7] + fvel*a[8]) + gsl_pow_4(musqr*fvel)*a[9])*legendrep(ell,muval));
  }

double doLPTmuint(void *params) {
  muint_params *paramsmu = (muint_params *) params;
  int ell = paramsmu->ellval;
  gsl_function TT;
  TT.function = &muint_integrand;
  TT.params = params;
  double i1, i1err;
  int mystat;
  size_t neval;
  mystat = gsl_integration_qag(&TT,0.,1.,epsabsI1,epsrelI1,limitI1,GSL_INTEG_GAUSS31,workspaceI1,&i1,&i1err);
  return(i1*(2.*((double) ell)+1.)); //this gets the normalization right.
  }

double QI1_integrand(double rval, void *params) {
  i1int_params *paramsi1 = (i1int_params *) params;
  double kval = paramsi1 ->kval;
  double xval = paramsi1 ->i2val;
  int qindx = paramsi1 ->n;
  return(pk(kval*rval)*pk(kval*sqrt(1.+rval*(rval-2.*xval)))*Qtilde(qindx,rval,xval));
  }

double QI1(void *params) {
  i1int_params *paramsi1 = (i1int_params *) params;
  int qindx = paramsi1 ->n;
  double xval = paramsi1 ->i2val;
  double kval = paramsi1 ->kval;
  gsl_function AA;
  AA.function = &QI1_integrand;
  AA.params = params;
  int I1intstatus;
  double i1, i1err;
  double i2, i2err;
  double i3, i3err;

  switch(qindx) {
    default:  //this works for 1, test for 2.
  //  case(1):
  //break into pieces around r=1, where Qtilde(r,x) peaks. 
  I1intstatus = gsl_integration_qag(&AA,0.,1.,epsabsI1,epsrelI1,limitI1,GSL_INTEG_GAUSS31,workspaceI1,&i1,&i1err);
  I1intstatus = gsl_integration_qag(&AA,1.,10.,epsabsI1,epsrelI1,limitI1,GSL_INTEG_GAUSS31,workspaceI1,&i2,&i2err);
  //to infty, don't need a lot of accuracy though.
  I1intstatus = gsl_integration_qagiu(&AA,10.,1.0,1.0e-2,limitI1,workspaceI1,&i3,&i3err);
  return(i1+i2+i3);
/*
    default:
  exit(1);
*/
  } //end switch on qindx.
  }

//integrate wrt y=kr instead of r.
double RI_integrand(double yval, void *params) {
  i2int_params *paramsi2 = (i2int_params *) params;
  double kval = paramsi2 -> kval;
  int rindx = paramsi2 -> n;
  if(kval < 1.0e-6) {
    kval = 1.0e-6; //avoid 0 in denom of Rtilde.
    }
  return(pk(yval)*Rtilde(rindx,yval/kval));
  }

double RI(double k, int n) {
  gsl_function CC;
  CC.function = &RI_integrand;
  i2int_params myp3 = {k, n};
  CC.params = &myp3;
  double intlimits[4];
  double i1,i1err;
  double i1tot = 0.;
  int I1intstatus;
  int i,j;
  intlimits[0] = 0.;
  intlimits[1] = 1.;
  intlimits[2] = 10.;
  intlimits[3] = k;
  for(i=0;i<=2;i++) {
    if(k < intlimits[i]) {
      for(j=2;j>=i;j--) {
        intlimits[j+1] = intlimits[j];
        }
      intlimits[i] = k;
      break;
      }
    }

  //take this out afterwards.
  for(i=0;i<=2;i++) {
    if(!(intlimits[i] <= intlimits[i+1])) {
      printf("tt %e %e %e\n",intlimits[i],intlimits[i+1],k);
      }
    assert(intlimits[i] <= intlimits[i+1]);
    }
  for(i=0;i<=2;i++) {
    I1intstatus = gsl_integration_qag(&CC,intlimits[i],intlimits[i+1],epsabsI2,epsrelI2,limitI2,GSL_INTEG_GAUSS31,workspaceI2,&i1,&i1err);
    i1tot += i1;
    }
  //to infty, don't need a lot of accuracy though.
  I1intstatus = gsl_integration_qagiu(&CC,intlimits[3],1.0,1.0e-2,limitI2,workspaceI2,&i1,&i1err);
  i1tot += i1;
  return(i1tot/k);  //since integral was over y=kr, need to multiply my result by 1/k.
  }

double QI2_integrand(double xval, void *params) {
  i2int_params *paramsi2 = (i2int_params *) params;
  int qindx = paramsi2 ->n;
  double kval = paramsi2 ->kval;

  i1int_params myp = {kval, xval, qindx};
  return(QI1(&myp)); 
  }

double QI2(double k, int n) {
  int I2intstatus;
  double i1, i1err;
  double *points;
  gsl_function BB;
  BB.function = &QI2_integrand;
  i2int_params myp2 = {k,n};
  BB.params = &myp2;
  if(n == 3 || n == 4 || n == 7 || n == 9 || n == 14 || n == 15) {
    //x = 1 is not well-behaved. at fixed x=1, there is a 1/(1-r)^2 divergence at r=1, which means the r integral won't converge.
    //for 7 and 9, divergence is 1/(1-r)
    points = (double *) malloc(2*sizeof(double));
    points[0] = -1.;
    points[1] = 1.;
    I2intstatus = gsl_integration_qagp(&BB,points,2,epsabsI2,epsrelI2,limitI2,workspaceI2,&i1,&i1err);
    return(i1);
    exit(1);
    }
  else {
    //  case(1):
    I2intstatus = gsl_integration_qag(&BB,-1.,1.,epsabsI2,epsrelI2,limitI2,GSL_INTEG_GAUSS15,workspaceI2,&i1,&i1err);
    return(i1);
    }
  exit(1);
  } //end QI2

double kNLsqr_integrand(double k, void *params) {
  return(pk(k));
  }

//we're defining sigNL^2 = 1./k_NL^2 = 1/6pi^2 * int dk P(k)
double getsigNLsqr() {
  if(ptsetup == 0) {
    ptintegralsetup();
    }
  gsl_set_error_handler_off();
  int I1intstatus;
  double i3,i3err;
  gsl_function DD;
  DD.function = &kNLsqr_integrand;
  DD.params = NULL;
  I1intstatus = gsl_integration_qagiu(&DD,0.,epsabsI2,epsrelI2,limitI1,workspaceI1,&i3,&i3err);
  return(i3/6./M_PI/M_PI);
  }

double Qkint(double k, int n) {
  assert(n >= 1 && n <= 17);  //last four are SPT integrals for Pdt and Ptt.
  if(ptsetup == 0) {
    ptintegralsetup();
    }
  gsl_set_error_handler_off();
  return QI2(k, n);
  }

double Rkint(double k, int n) {
  assert(n >= 1 && n <= 8); //last two are SPT integrals for Pdt and Ptt. 
  if(ptsetup == 0) {
    ptintegralsetup();
    }
  gsl_set_error_handler_off();
  return RI(k,n);
  }

//use fine grid on low k (where BAO is, etc), and coarser grid at high k.

int setupSPTvelsplines() {
  int j;
  assert(lenkvelPT > 0);
  splinevinlist = (gsl_spline **) malloc((2)*sizeof(gsl_spline *));
  accvinlist = (gsl_interp_accel **) malloc((2)*sizeof(gsl_interp_accel *));
  splinesig2list = (gsl_spline **) malloc((5)*sizeof(gsl_spline *));
  accsig2list = (gsl_interp_accel **) malloc((5)*sizeof(gsl_interp_accel *));
  for(j=0;j<2;j++) {
    splinevinlist[j] = gsl_spline_alloc(gsl_interp_cspline,lenkvelPT);
    accvinlist[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinevinlist[j], klistvelPT, vinlist[j], lenkvelPT);    
    }
  for(j=0;j<5;j++) {
    splinesig2list[j] = gsl_spline_alloc(gsl_interp_cspline,lenkvelPT);
    accsig2list[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinesig2list[j], klistvelPT, sig2list[j], lenkvelPT);    
    }
  return 0;
  }

int doSPTvelints(double kmin, double kmax1, double kmax2, double kmax3, double dk1, double dk2, double dk3) {
  time_t time1, time2;
  (void) time(&time1);
  if(dk1 < 1.0e-5*0.7) {
    dk1 = 1.0e-3*0.7;
    }
  if(dk2 < 1.0e-5*0.7) {
    dk2 = 1.0e-3*0.7;
    }
  if(dk3 < 1.0e-5*0.7) {
    dk3 = 1.0e-3*0.7;
    }
  int l1, l2, l3;
  l1 = ((int) (floor((kmax1 - kmin)/dk1)+1.));
  //weird bug -- got a crazy value at k=7.1652, but not if i go around it.
  //for now just stop at lower k, but be aware that the integator can do weird things at high k.
  //l2 = ((int) (floor((kmax2 - (kmax1+dk2))/dk2)+2.)); //+2 just to be extra safe.
  //l3 = ((int) (floor((kmax3 - (kmax2+dk3))/dk3)+2.)); //+2 just to be extra safe.
  l2 = ((int) (floor((kmax2 - (kmax1+dk2))/dk2)+1.));
  l3 = ((int) (floor((kmax3 - (kmax2+dk3))/dk3)+1.));
  if(kmax2 < kmax1) {
    l2 = 0;
    }
  if(kmax3 < kmax2) {
    l3 = 0;
    }
  lenkvelPT = l1 + l2 + l3;
  klistvelPT = (double *) malloc(sizeof(double)*lenkvelPT);
  //logklistvelPT = (double *) malloc(sizeof(double)*lenkvelPT);
  vinlist = malloc2ddouble(2,lenkvelPT); //malloc 2d array.
  sig2list = malloc2ddouble(5,lenkvelPT);
  double k, kscalefac;
  int i;
  for(i=0;i<lenkvelPT;i++) {
    #ifdef VERBOSE
    if(i%10==0) {
      (void) time(&time2);
      fprintf(stderr,"doSPTints %d of %d, delta t = %.3e\n",i,lenkvelPT,difftime(time2,time1));
      fflush(stderr);
      }
    #endif
    if(i == 0) {
      k = kmin;
      }
    else {
      if(klistvelPT[i-1] + dk1 > kmax1) {
        if(klistvelPT[i-1]+dk2 > kmax2) {
          k = klistvelPT[i-1] + dk3;
          }
        else {
          k = klistvelPT[i-1] + dk2;
          }
        }
      else {
        k = klistvelPT[i-1] + dk1;
        }
      }
    klistvelPT[i] = k;
    kscalefac = gsl_pow_3(k)/4./M_PI/M_PI;
    vinlist[0][i] = ((Qkint(k,14) + Rkint(k,3)*pk(k))*kscalefac)*k;  //extra k from derivative in vin.`
    sig2list[0][i] = ((Qkint(k,15) + Rkint(k,4)*pk(k))*kscalefac);
    vinlist[1][i] = (Qkint(k,16) + Rkint(k,5)*pk(k))*kscalefac*k; //still needs to be multiplied by -2b^2f/pi^2 and r integral done.  velocity has an extra k from the derivative.
    sig2list[1][i] = Qkint(k,17)*gsl_pow_3(k)/M_PI/M_PI; //still needs to be multiplied by f^2/pi^2*b  Eqn A14.
    sig2list[2][i] = Rkint(k,6)*pk(k)*gsl_pow_3(k)/M_PI/M_PI; //still needs to be multiplied by f^2/pi^2*b  Eqn A12 top.
    sig2list[3][i] = Rkint(k,7)*pk(k)*gsl_pow_3(k)/M_PI/M_PI; //still needs to be multiplied by f^2/pi^2*b  Eqn A12 bottom
    sig2list[4][i] = Rkint(k,8)*pk(k)*gsl_pow_3(k)/M_PI/M_PI; //still needs to be multiplied by f^2/pi^2*b  Eqn A15.  switched sign convention from iset1.c on this integral, careful!
    }
  setupSPTvelsplines();
  return(0);
  }

int setupQRsplines() {
  int j;
  splineqlist = (gsl_spline **) malloc((13)*sizeof(gsl_spline *));
  accqlist = (gsl_interp_accel **) malloc((13)*sizeof(gsl_interp_accel *));
  splinerlist = (gsl_spline **) malloc((2)*sizeof(gsl_spline *));
  accrlist = (gsl_interp_accel **) malloc((2)*sizeof(gsl_interp_accel *));
  for(j=0;j<13;j++) {
    splineqlist[j] = gsl_spline_alloc(gsl_interp_cspline,lenkqr);
    accqlist[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splineqlist[j], klistqr, qlist[j], lenkqr);    
    }
  for(j=0;j<2;j++) {
    splinerlist[j] = gsl_spline_alloc(gsl_interp_cspline,lenkqr);
    accrlist[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinerlist[j], klistqr, rlist[j], lenkqr);    
    }
  return 0;
  }


double interpvinlist(double k, int jj) {
  if(k < klistvelPT[0]) {
    //linearly interpolate to 0.
    return(vinlist[jj][0]*k/klistvelPT[0]);
    }
  if(k > klistvelPT[lenkvelPT-1]) {  //make sure we evaluate out to all k values we care about.
    assert(k < 1.05*klistvelPT[lenkvelPT-1]);
    // use linear interpolation outside cspline boundaries.
    double slope = (vinlist[jj][lenkvelPT-1] - vinlist[jj][lenkvelPT-2])/(klistvelPT[lenkvelPT-1] - klistvelPT[lenkvelPT - 2]);
    return (vinlist[jj][lenkvelPT-1] + slope * (k - klistvelPT[lenkvelPT-1]));
    }
  return(gsl_spline_eval(splinevinlist[jj],k,accvinlist[jj]));
  }

double interpsig2list(double k, int jj) {
  if(k < klistvelPT[0]) {
    //linearly interpolate to 0.
    return(sig2list[jj][0]*k/klistvelPT[0]);
    }
  if(k > klistvelPT[lenkvelPT-1]) {  //make sure we evaluate out to all k values we care about.
    assert(k < 1.05*klistvelPT[lenkvelPT-1]);
    // use linear interpolation outside cspline boundaries.
    double slope = (sig2list[jj][lenkvelPT-1] - sig2list[jj][lenkvelPT-2])/(klistvelPT[lenkvelPT-1] - klistvelPT[lenkvelPT - 2]);
    return (sig2list[jj][lenkvelPT-1] + slope * (k - klistvelPT[lenkvelPT-1]));
    }
  return(gsl_spline_eval(splinesig2list[jj],k,accsig2list[jj]));
  }

double interpQ(double k, int jj) {
  if(k < klistqr[0]) {
    //linearly interpolate to 0.
    return(qlist[jj][0]*k/klistqr[0]);
    }
  if(k > klistqr[lenkqr-1]) {  //make sure we evaluate out to all k values we care about.
    assert(k < 1.01*klistqr[lenkqr-1]);
    return 0;
    }
  return(gsl_spline_eval(splineqlist[jj],k,accqlist[jj]));
  }

double interpR(double k, int jj) {
  if(k < klistqr[0]) {
    //linearly interpolate to 0.
    return(rlist[jj][0]*k/klistqr[0]);
    }
  if(k > klistqr[lenkqr-1]) {  //make sure we evaluate out to all k values we care about.
    assert(k < 1.01*klistqr[lenkqr-1]);
    return 0;
    }
  return(gsl_spline_eval(splinerlist[jj],k,accrlist[jj]));
  }

//get Q integrals for every k in the input list.
int doQRints(double kmin, double kmaxqr1, double kmaxqr2, double dkqr1, double dkqr2) {
  time_t time1, time2;
  (void) time(&time1);
  if(dkqr1 < 1.0e-5*0.7) {
    dkqr1 = 1.0e-3*0.7;
    }
  if(dkqr2 < 1.0e-5*0.7) {
    dkqr2 = 1.0e-3*0.7;
    }
  int l1, l2;
  l1 = ((int) (floor((kmaxqr1 - kmin)/dkqr1)+1.));
  l2 = ((int) (floor((kmaxqr2 - (kmaxqr1+dkqr2))/dkqr2)+2.)); //+2 just to be extra safe.
  if(kmaxqr2 < kmaxqr1) {
    l2 = 0;
    }
  lenkqr = l1 + l2;
  klistqr = (double *) malloc(sizeof(double)*lenkqr);
  qlist = malloc2ddouble(13,lenkqr); //malloc 2d array.
  rlist = malloc2ddouble(2,lenkqr); //malloc 2d array.
  double k;
  int i,j;
  double kscalefac;
  for(i=0;i<lenkqr;i++) {
    #ifdef VERBOSE
    if(i%10==0) {
      (void) time(&time2);
      fprintf(stderr,"doQRints %d of %d, delta t = %.3e\n",i,lenkqr,difftime(time2,time1));
      fflush(stderr);
      }
    #endif
    if(i == 0) {
      k = kmin;
      }
    else {
      if(klistqr[i-1] + dkqr1 > kmaxqr1) {
        k = klistqr[i-1] + dkqr2;
        }
      else {
        k = klistqr[i-1] + dkqr1;
        }
      }
    klistqr[i] = k;
    kscalefac = gsl_pow_3(k)/4./M_PI/M_PI;
    for(j=0;j<13;j++) {
      qlist[j][i]=Qkint(k,j+1)*kscalefac; //now qlist is exactly Eqn A39
      }
    for(j=0;j<2;j++) {
      rlist[j][i]=Rkint(k,j+1)*kscalefac*pk(k); //now rlist is exactly Eqn A40
      }
    }
  setupQRsplines();
  return(0);
  }


double vinSPTkint(double k) {
  if(k < klistvel[0]) {
    assert(b1Lhold > -50);
    return(pk(k)*k*(1.+b1Lhold));
    }
  if(k > klistvel[lenkvel-1]) {
    return  (exp(logklistvel[lenkvel-1] + nasympvISPT*(log(k) - logklistvel[lenkvel-1])));
    }
  return (exp(gsl_spline_eval(splinevISPT,log(k),accvISPT)));
  }

double sigsqrSPTkint(double k, int whichint) {
  if(k < klistvel[0]) {
    if(whichint == 0) {
      return(pk(k));
      }
    else {
      return 0.;
      }
    }
  if(k > klistvel[lenkvel-1]) {
    if(whichint == 0) {
      return  (exp(logklistvel[lenkvel-1] + nasympsigsqrSPT0*(log(k) - logklistvel[lenkvel-1])));
      }
    if(whichint == 1) {
      return  (exp(logklistvel[lenkvel-1] + nasympsigsqrSPT1*(log(k) - logklistvel[lenkvel-1])));
      }
    if(whichint == 2) {
      return  (exp(logklistvel[lenkvel-1] + nasympsigsqrSPT2*(log(k) - logklistvel[lenkvel-1])));
      }
    exit(1);
    }
  if(whichint == 0) {
    return (exp(gsl_spline_eval(splinesigsqrSPT0,log(k),accsigsqrSPT0)));
    }
  if(whichint == 1) {
    return (exp(gsl_spline_eval(splinesigsqrSPT1,log(k),accsigsqrSPT1)));
    }
  if(whichint == 2) {
    return (exp(gsl_spline_eval(splinesigsqrSPT2,log(k),accsigsqrSPT2)));
    }
  printf("sigsqrSPTkint out of range\n");
  exit(1);
  }

int getkintvsigrealSPT(double b1L, double b2L, double kmin, double kmax1, double kmax2, double dk1, double dk2) {
  if(b1Lhold > 0) {
    assert(fabs(b1Lhold - b1L) < 1.0e-6);
    }
  else {
    b1Lhold = b1L;
    }

//assign klistvel here, more finely spaced than PT integrals.
  if(dk1 < 1.0e-5*0.7) {
    dk1 = 1.0e-3*0.7;
    }
  if(dk2 < 1.0e-5*0.7) {
    dk2 = 1.0e-3*0.7;
    }
  int l1, l2;
  l1 = ((int) (floor((kmax1 - kmin)/dk1)+1.));
  l2 = ((int) (floor((kmax2 - (kmax1+dk2))/dk2)+2.)); //+2 just to be extra safe.
  if(kmax2 < kmax1) {
    l2 = 0;
    }
  lenkvel = l1 + l2;
  klistvel = (double *) malloc(sizeof(double)*lenkvel);
  logklistvel = (double *) malloc(sizeof(double)*lenkvel);
  double k;
  int i;
  for(i=0;i<lenkvel;i++) {
    if(i == 0) {
      k = kmin;
      }
    else {
      if(klistvel[i-1] + dk1 > kmax1) {
        k = klistvel[i-1] + dk2;
        }
      else {
        k = klistvel[i-1] + dk1;
        }
      }
    klistvel[i] = k;
    logklistvel[i] = log(k);
    }
  logvISPT = (double *) malloc(sizeof(double)*lenkvel);
  logsigsqrSPT0 = (double *) malloc(sizeof(double)*lenkvel);
  logsigsqrSPT1 = (double *) malloc(sizeof(double)*lenkvel);
  logsigsqrSPT2 = (double *) malloc(sizeof(double)*lenkvel);
  for(i=0;i<lenkvel;i++) {
    k = klistvel[i];
    logvISPT[i] = (1.+b1L)*(k*pk(k)+interpvinlist(k,0)) + 2.*(1.+b1L)*(1.+b1L)*interpvinlist(k,1);
    logsigsqrSPT0[i] = pk(k) + interpsig2list(k,0); 
    logsigsqrSPT1[i] = (1.+b1L)*(interpsig2list(k,1) + interpsig2list(k,3) + interpsig2list(k,4));
    logsigsqrSPT2[i] = (1.+b1L)*(interpsig2list(k,2));

    //these interpolations are csplines and behave wacky near the boundary.
    if((logvISPT[i] <= 0.) || isnan(logvISPT[i])) {
      assert(0==1); //are we ever here?  i hope changes to interpolations mean we're not.
      assert(i>0);
      logvISPT[i] = logvISPT[i-1];
      }
    if((logsigsqrSPT2[i] <= 0.) || isnan(logsigsqrSPT2[i])) {
      assert(0==1); //are we ever here?  i hope changes to interpolations mean we're not.
      assert(i>0);
      logsigsqrSPT2[i] = logsigsqrSPT2[i-1];
      }
    if((logsigsqrSPT1[i] <= 0.) || isnan(logsigsqrSPT1[i])) {
      assert(0==1); //are we ever here?  i hope changes to interpolations mean we're not.
      assert(i>0);
      logsigsqrSPT1[i] = logsigsqrSPT1[i-1];
      }
    if((logsigsqrSPT0[i] <= 0.) || isnan(logsigsqrSPT0[i])) {
      assert(0==1); //are we ever here?  i hope changes to interpolations mean we're not.
      assert(i>0);
      logsigsqrSPT0[i] = logsigsqrSPT0[i-1];
      }

    assert(logvISPT[i] > 0.);
    assert(logsigsqrSPT0[i] > 0.);
    assert(logsigsqrSPT1[i] > 0.);
    assert(logsigsqrSPT2[i] > 0.);
    logvISPT[i] = log(logvISPT[i]);
    logsigsqrSPT0[i] = log(logsigsqrSPT0[i]);
    logsigsqrSPT1[i] = log(logsigsqrSPT1[i]);
    logsigsqrSPT2[i] = log(logsigsqrSPT2[i]);
    }

  accvISPT = gsl_interp_accel_alloc ();
  splinevISPT = gsl_spline_alloc (gsl_interp_linear, lenkvel); 
  gsl_spline_init (splinevISPT, logklistvel, logvISPT, lenkvel);
  nasympvISPT = gsl_spline_eval_deriv(splinevISPT, logklistvel[lenkvel-1], accvISPT);
  #ifdef VERBOSE
  printf("vISPT asymp: %e\n",nasympvISPT);
  #endif

  accsigsqrSPT0 = gsl_interp_accel_alloc ();
  splinesigsqrSPT0 = gsl_spline_alloc (gsl_interp_linear, lenkvel); 
  gsl_spline_init (splinesigsqrSPT0, logklistvel, logsigsqrSPT0, lenkvel);
  nasympsigsqrSPT0 = gsl_spline_eval_deriv(splinesigsqrSPT0, logklistvel[lenkvel-1], accsigsqrSPT0);
  #ifdef VERBOSE
  printf("sigsqrSPT0 asymp: %e\n",nasympsigsqrSPT0);
  #endif

  accsigsqrSPT1 = gsl_interp_accel_alloc ();
  splinesigsqrSPT1 = gsl_spline_alloc (gsl_interp_linear, lenkvel); 
  gsl_spline_init (splinesigsqrSPT1, logklistvel, logsigsqrSPT1, lenkvel);
  nasympsigsqrSPT1 = gsl_spline_eval_deriv(splinesigsqrSPT1, logklistvel[lenkvel-1], accsigsqrSPT1);
  #ifdef VERBOSE
  printf("sigsqrSPT1 asymp: %e\n",nasympsigsqrSPT1);
  #endif

  accsigsqrSPT2 = gsl_interp_accel_alloc ();
  splinesigsqrSPT2 = gsl_spline_alloc (gsl_interp_linear, lenkvel); 
  gsl_spline_init (splinesigsqrSPT2, logklistvel, logsigsqrSPT2, lenkvel);
  nasympsigsqrSPT2 = gsl_spline_eval_deriv(splinesigsqrSPT2, logklistvel[lenkvel-1], accsigsqrSPT2);
  #ifdef VERBOSE
  printf("sigsqrSPT2 asymp: %e\n",nasympsigsqrSPT2);
  #endif

  //sanity check.
  for(i=0;i<lenkvel;i++) {
    assert(fabs(exp(logvISPT[i]) - vinSPTkint(klistvel[i]))/exp(logvISPT[i]) < 1.0e-5);
    assert(fabs(exp(logsigsqrSPT0[i]) - sigsqrSPTkint(klistvel[i],0))/exp(logsigsqrSPT0[i]) < 1.0e-5);
    assert(fabs(exp(logsigsqrSPT1[i]) - sigsqrSPTkint(klistvel[i],1))/exp(logsigsqrSPT1[i]) < 1.0e-5);
    assert(fabs(exp(logsigsqrSPT2[i]) - sigsqrSPTkint(klistvel[i],2))/exp(logsigsqrSPT2[i]) < 1.0e-5);
    }
  #ifdef VERBOSE
  printf("vinSPTkint and sigsqr spline passed sanity test.\n");
  #endif
  return 0;
  }

int getpkrealLPT(double b1L, double b2L, double kminLPT, double kmaxLPT, double dkLPT) {
  int i;
  if(lenk == 0) {
    lenk = ((int) (floor((kmaxLPT - kminLPT)/dkLPT)+1.));
    klist = (double *) malloc(sizeof(double)*lenk);
    for(i=0;i<lenk;i++) {
      klist[i] = kminLPT + ((double) i)*dkLPT;
      }
    }
  
  pkLPTreal = (double *) malloc(sizeof(double)*lenk);
  double sigNLsqr = getsigNLsqr();
  #ifdef VERBOSE
  printf("sigNL, kNL = %e %e\n",sqrt(sigNLsqr),1./sqrt(sigNLsqr));
  #endif
  double e00;
  double prefac,k;
  for(i=0;i<lenk;i++) {
    prefac = exp(-klist[i]*klist[i]*sigNLsqr);
    k = klist[i];
    //remember Q_i in Matsubara is saved here as qlist[i-1]

    e00 = (9./98.*interpQ(k,0) + 3./7.*interpQ(k,1) + 0.5*interpQ(k,2) + 10./21.*interpR(k,0) + 6./7.*interpR(k,1))
    + b1L*(6./7.*interpQ(k,4) + 2.*interpQ(k,6) + 4./3.*interpR(k,0) + 12./7.*interpR(k,1))
    + b2L*(3./7.*interpQ(k,7) + interpQ(k,8))
    + b1L*b1L*(interpQ(k,8) + interpQ(k,10) + 6./7.*(interpR(k,0) + interpR(k,1)))
    + 2.*b1L*b2L*(interpQ(k,11))
    + 0.5*b2L*b2L*(interpQ(k,12));
    pkLPTreal[i] = prefac*((1.+b1L)*(1.+b1L)*pk(k)+e00);
    }

//copied from linear pk spline.
  accpkLPTreal = gsl_interp_accel_alloc ();
  splinepkLPTreal = gsl_spline_alloc (gsl_interp_linear, lenk); 
  gsl_spline_init (splinepkLPTreal, klist, pkLPTreal, lenk);
  return 0;
  }

int getpkredLPT(double b1L, double b2L, double fvel, int ellmax, double kminLPT, double kmaxLPT, double dkLPT) {

  int i,j;
  if(b1Lhold > 0) {
    assert(fabs(b1Lhold - b1L) < 1.0e-6);
    }
  else {
    b1Lhold = b1L;
    }
  if(fvelhold > 0) {
    assert(fabs(fvelhold - fvel) < 1.0e-6);
    }
  else {
    fvelhold = fvel;
    }
  if(lenk == 0) {
    lenk = ((int) (floor((kmaxLPT - kminLPT)/dkLPT)+1.));
    klist = (double *) malloc(sizeof(double)*lenk);
    for(i=0;i<lenk;i++) {
      klist[i] = kminLPT + ((double) i)*dkLPT;
      }
    }


  pkLPTredellmax = ellmax;
  //opposite order of rest of 2d arrays; but we want to be able to interpolate over k at each ell.
  pkLPTred = malloc2ddouble(ellmax/2+1,lenk);
  double sigNLsqr = getsigNLsqr();
  #ifdef VERBOSE
  printf("sigNL, kNL = %e %e\n",sqrt(sigNLsqr),1./sqrt(sigNLsqr));
  #endif
  double e00, e11, e12, e22, e23, e24, e33, e34, e44;
  double prefac,k;
  //could probably integrate analytically, but let's do it numerically for now.
  double aa[10];
  aa[0] = sigNLsqr;

  for(i=0;i<lenk;i++) {
    k = klist[i];
    e00 = (9./98.*interpQ(k,0) + 3./7.*interpQ(k,1) + 0.5*interpQ(k,2) + 10./21.*interpR(k,0) + 6./7.*interpR(k,1))
    + b1L*(6./7.*interpQ(k,4) + 2.*interpQ(k,6) + 4./3.*interpR(k,0) + 12./7.*interpR(k,1))
    + b2L*(3./7.*interpQ(k,7) + interpQ(k,8))
    + b1L*b1L*(interpQ(k,8) + interpQ(k,10) + 6./7.*(interpR(k,0) + interpR(k,1)))
    + 2.*b1L*b2L*(interpQ(k,11))
    + 0.5*b2L*b2L*(interpQ(k,12));

    e11 = (18./49.*interpQ(k,0) + 12./7.*interpQ(k,1) + 2.*interpQ(k,2) + 40./21.*interpR(k,0) + 24./7.*interpR(k,1))
    + b1L*(18./7.*interpQ(k,4) + 6.*interpQ(k,6) + 4.*interpR(k,0) + 36./7.*interpR(k,1)) 
    + b2L*(6./7.*interpQ(k,7) + 2.*interpQ(k,8))
    + b1L*b1L*(2.*interpQ(k,8) + 2.*interpQ(k,10) + 12./7.*(interpR(k,0) + interpR(k,1)))
    + 2.*b1L*b2L*(interpQ(k,11));

    e12 = (-3./14.*interpQ(k,0) - 1.5*interpQ(k,1) + 0.25*interpQ(k,3) - 6./7.*interpR(k,0))
    + b1L*(interpQ(k,5) - 6./7.*interpR(k,0))
    - 0.5*b2L*interpQ(k,7)
    - 0.5*b1L*b1L*(interpQ(k,7) - interpQ(k,9));



    e22 = (57./98.*interpQ(k,0) + 51./14.*interpQ(k,1) + 3.*interpQ(k,2) -0.25*interpQ(k,3) + 16./7.*interpR(k,0) + 30./7.*interpR(k,1))
    + b1L*(12./7.*interpQ(k,4) - interpQ(k,5) + 6.*interpQ(k,6) + 18./7.*interpR(k,0) + 24./7.*interpR(k,1))
    + b2L*(0.5*interpQ(k,7) + interpQ(k,8))
    + b1L*b1L*(0.5*interpQ(k,7) + interpQ(k,8) -0.5*interpQ(k,9) + interpQ(k,10));


    e23 = -3./7.*interpQ(k,0) -3.*interpQ(k,1) + 0.5*interpQ(k,3) - 6./7.*interpR(k,0) + b1L*interpQ(k,5);

    e24 = 3./16.*interpQ(k,0);

    e33 = 3./7.*interpQ(k,0) + 27./7.*interpQ(k,1) + 2.*interpQ(k,2) - 0.5*interpQ(k,3) + 6./7.*interpR(k,0) + 12./7.*interpR(k,1) 
    + b1L*(-interpQ(k,5) + 2.*interpQ(k,6));

    e34 = -3./8.*interpQ(k,0)-1.5*interpQ(k,1)+0.25*interpQ(k,3);

    e44 = 3./16.*interpQ(k,0) + 1.5*interpQ(k,1) + 0.5*interpQ(k,2) - 0.25*interpQ(k,3);


    aa[1] = e00;
    aa[2] = e11;
    aa[3] = e12;
    aa[4] = e22;
    aa[5] = e23;
    aa[6] = e24;
    aa[7] = e33;
    aa[8] = e34;
    aa[9] = e44;

    for(j=0;j<=ellmax/2;j++) {
      muint_params mypmu = {aa, klist[i], b1L, fvel, j*2};
      pkLPTred[j][i] = doLPTmuint(&mypmu);
      }
    }

  //now set up splines.

  splinepkLPTred = (gsl_spline **) malloc((ellmax/2+1)*sizeof(gsl_spline *));
  accpkLPTred = (gsl_interp_accel **) malloc((ellmax/2+1)*sizeof(gsl_interp_accel *));
  for(j=0;j<=ellmax/2;j++) {
    splinepkLPTred[j] = gsl_spline_alloc(gsl_interp_linear,lenk);
    accpkLPTred[j] = gsl_interp_accel_alloc();
    gsl_spline_init (splinepkLPTred[j], klist, pkLPTred[j], lenk);    
    }
  return 0;
  }

int printpkredLPT(char *pkoutfname) {
  int i,j;
  FILE *ofp;
  ofp = open_file_write(pkoutfname);
  for(i=0;i<lenk;i++) {
    fprintf(ofp,"%e ",klist[i]);
    for(j=0;j<=pkLPTredellmax/2;j++) {
      fprintf(ofp,"%e ",pkLPTred[j][i]);
      assert(fabs(pkLPTred[j][i]-pkredLPT(klist[i],2*j))/fabs(max(pkLPTred[j][i],1.0)) < 1.0e-6);
      }
    fprintf(ofp,"\n");
    }
  fclose(ofp);
  return 0;
  }

double pkLPT(double k) {
  if(k < klist[0]) {
    //use linear.
    return(pk(k));
    }
  if(k > klist[lenk-1]) {
    return 0.;
    }
  return(gsl_spline_eval(splinepkLPTreal,k,accpkLPTreal));
  }


double pkredLPT(double k, int ell) {
  assert(ell <= pkLPTredellmax);
  if(k < klist[0]) {
    //use linear theory.
    if(ell == 0) {
      return((gsl_pow_2(1.+b1Lhold)+2./3.*(1.+b1Lhold)*fvelhold+1./5.*fvelhold*fvelhold)*pk(k));
      }
    if(ell == 2) {
      return((4./3.*(1.+b1Lhold)*fvelhold+4./7.*fvelhold*fvelhold)*pk(k));
      }
    if(ell == 4) {
      return((8./35.*fvelhold*fvelhold)*pk(k));
      }
    return 0.; //all other multipoles 0 in linear theory.
    }
  if(k > klist[lenk-1]) {
    return 0.;
    }
  return(gsl_spline_eval(splinepkLPTred[ell/2],k,accpkLPTred[ell/2]));
  }


int printpkreal(char *pkoutfname) {
  int i,j;
  FILE *ofp;
  ofp = open_file_write(pkoutfname);
  for(i=0;i<lenk;i++) {
    fprintf(ofp,"%.4e %.4e %.4e\n",klist[i],pk(klist[i]),pkLPTreal[i]);
    }
  fclose(ofp);
  }

int printQRints(char *QRfname) {
  int i,j;
  FILE *ofp;
  ofp = open_file_write(QRfname);
  for(i=0;i<lenkqr;i++) {
    fprintf(ofp,"%e ",klistqr[i]);
    for(j=0;j<13;j++) {
      fprintf(ofp,"%e ",qlist[j][i]);
      }
    for(j=0;j<2;j++) {
      fprintf(ofp,"%e ",rlist[j][i]);
      }
    fprintf(ofp,"\n");
    }
  fclose(ofp);
  return 0;
  }

int printSPTvelints(char *SPTfname) {
  int i,j;
  FILE *ofp;
  ofp = open_file_write(SPTfname);
  for(i=0;i<lenkvelPT;i++) {
    fprintf(ofp,"%e ",klistvelPT[i]);
    for(j=0;j<2;j++) {
      fprintf(ofp,"%e ",vinlist[j][i]);
      }
    for(j=0;j<5;j++) {
      fprintf(ofp,"%e ",sig2list[j][i]);
      } 
    fprintf(ofp,"\n");
    }
  fclose(ofp);
  return 0;
  }

int readSPTvelints(char *SPTfname) {
  FILE *ifpspt;
  char line[MAXLINELEN];
  int headercount,linecount;
  ifpspt = open_file_read(SPTfname);
  get_file_length(ifpspt,&headercount,&linecount);
  assert(headercount == 0); //may want to add a header later.
  lenkvelPT = linecount - headercount;
  klistvelPT = (double *) malloc(sizeof(double)*lenkvelPT);
  //logklistvel = (double *) malloc(sizeof(double)*lenkvelPT);
  vinlist = malloc2ddouble(2,lenkvelPT); //malloc 2d array.
  sig2list = malloc2ddouble(5,lenkvelPT);

  int i,j;
  //skip over header.
  for(i=0;i<headercount;i++) {
    fgets(line,MAXLINELEN,ifpspt);
    }

  //read in data.
  for(i=0;i<lenkvelPT;i++) {
    fscanf(ifpspt,"%le",&(klistvelPT[i]));
    for(j=0;j<2;j++) {
      fscanf(ifpspt,"%le",&(vinlist[j][i]));
      }
    for(j=0;j<5;j++) {
      fscanf(ifpspt,"%le",&(sig2list[j][i]));
      }
    }
  fclose(ifpspt);
  setupSPTvelsplines();
  return(0);
  }

int readQRints(char *QRfname) {
  FILE *ifpqr;
  char line[MAXLINELEN];
  int headercount,linecount;
  ifpqr = open_file_read(QRfname);
  get_file_length(ifpqr,&headercount,&linecount);
  assert(headercount == 0); //may want to add a header later.
  lenkqr = linecount - headercount;
  klistqr = (double *) malloc(sizeof(double)*lenkqr);
  qlist = malloc2ddouble(13,lenkqr); //malloc 2d array.
  rlist = malloc2ddouble(2,lenkqr); //malloc 2d array.
  int i,j;
  //skip over header.
  for(i=0;i<headercount;i++) {
    fgets(line,MAXLINELEN,ifpqr);
    }

  //read in data.
  for(i=0;i<lenkqr;i++) {
    fscanf(ifpqr,"%le ",&(klistqr[i]));
    for(j=0;j<13;j++) {
      fscanf(ifpqr,"%le ",&(qlist[j][i]));
      }
    for(j=0;j<2;j++) {
      fscanf(ifpqr,"%le ",&(rlist[j][i]));
      }
    fscanf(ifpqr,"\n");
    }
  fclose(ifpqr);
  setupQRsplines();
  return(0);
  }

int freeptraw() {
  ptintegralfree();
  int j; 

  if(klistqr) {
    free(klistqr);
    }
  if(qlist) {
    malloc2dfree(qlist);
    }
  if(rlist) {
    malloc2dfree(rlist);
    }
  if(klist) {
    free(klist);
    }
  if(pkLPTreal) {
    free(pkLPTreal);
    }
  if(pkLPTred) {
    malloc2dfree(pkLPTred);
    }
  if(splinepkLPTred[0]) {
    for(j=0;j<=pkLPTredellmax/2;j++) {
      gsl_spline_free(splinepkLPTred[j]);
      gsl_interp_accel_free(accpkLPTred[j]);
      }
    free(splinepkLPTred);
    free(accpkLPTred); 
    }

  if(splineqlist[0]) {
    for(j=0;j<13;j++) {
      gsl_spline_free(splineqlist[j]);
      gsl_interp_accel_free(accqlist[j]);
      }
    free(splineqlist);
    free(accqlist); 
    }
  if(splinerlist[0]) {
    for(j=0;j<2;j++) {
      gsl_spline_free(splinerlist[j]);
      gsl_interp_accel_free(accrlist[j]);
      }
    free(splinerlist);
    free(accrlist); 
    }

  if(klistvel) {
    free(klistvel);
    }
  if(klistvelPT) {
    free(klistvelPT);
    }
  if(vinlist) {
    malloc2dfree(vinlist);
    }
  if(sig2list) {
    malloc2dfree(sig2list);
    }
  if(splinevinlist[0]) {
    for(j=0;j<2;j++) {
      gsl_spline_free(splinevinlist[j]);
      gsl_interp_accel_free(accvinlist[j]);
      }
    free(splinevinlist);
    free(accvinlist); 
    }
  if(splinesig2list[0]) {
    for(j=0;j<5;j++) {
      gsl_spline_free(splinesig2list[j]);
      gsl_interp_accel_free(accsig2list[j]);
      }
    free(splinesig2list);
    free(accsig2list); 
    }

  //free splines.
  if(splinepkLPTreal) {
    gsl_spline_free(splinepkLPTreal);
    gsl_interp_accel_free(accpkLPTreal);
    }
  }
