#include "linearpk.h"

//data structures holding linear power spectrum.
  static double *lnkin, *lnpkin;
  static int npk;
  static double nasymplow, nasymphigh;

 //power spectrum spline.
  static gsl_interp_accel *accpk;
  static gsl_spline *splinepk;

/***********************
  NOTES

  1) may want to play with spline type (currently set to linear in log-log space).

***********************/

int getpkspline(char *pkfname, double pkhval, double pkampscale) {
  FILE *ifppk;
  char line[MAXLINELEN];
  int i;
  int headercount,linecount;
  ifppk = open_file_read(pkfname);
  get_file_length(ifppk,&headercount,&linecount);
  npk = linecount - headercount;
  double kin, pkin;
  lnkin = (double *) malloc(npk*sizeof(double));
  lnpkin = (double *) malloc(npk*sizeof(double));

  //skip over header.
  for(i=0;i<headercount;i++) {
    fgets(line,MAXLINELEN,ifppk);
    }
  //read in pk data.
  for(i=0;i<npk;i++) {
    fscanf(ifppk,"%le %le\n",&(kin),&(pkin));
    lnkin[i] = log(kin*pkhval);
    lnpkin[i] = log(pkin*gsl_pow_3(1./pkhval)*pkampscale);
    }

  accpk = gsl_interp_accel_alloc ();
  splinepk = gsl_spline_alloc (gsl_interp_cspline, npk); 
  gsl_spline_init (splinepk, lnkin, lnpkin, npk);
  nasymphigh = gsl_spline_eval_deriv(splinepk, lnkin[npk-1], accpk);
  nasymplow = gsl_spline_eval_deriv(splinepk, lnkin[0], accpk);
  #ifdef VERBOSE
  printf("nasymps: %e %e\n",nasymphigh,nasymplow);
  #endif
  //sanity checks on asymptotic n values.
  assert(fabs(1-nasymplow) < 0.5);
  assert(nasymphigh < -2.);  //make sure we get out to reasonably high k.
  return 0;
  }

//use spline in region where defined; use power law extrapolation elsewhere.
double pk(double kval) {  //!!! k in units Mpc^-1.
  double lnkval = log(kval);
  if(lnkval >= lnkin[0] && lnkval <= lnkin[npk-1]) {
    return (exp(gsl_spline_eval(splinepk,lnkval,accpk)));
    }
  else {
    if(lnkval < lnkin[0]) {
      return (exp(lnpkin[0] + nasymplow*(lnkval - lnkin[0])));
      }
    else {
      return (exp(lnpkin[npk-1] + nasymphigh*(lnkval - lnkin[npk-1])));
      }
    }
  }

int pkfree() {
  gsl_spline_free(splinepk);
  gsl_interp_accel_free(accpk);
  free(lnkin);
  free(lnpkin);
  return 0;
  }


