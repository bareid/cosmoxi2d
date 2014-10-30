#include "besselroots.h"

int static besselsetup = 0;
static double besselroots[BMAX][2];

int static besselsetupA13 = 0;
static double besselrootsA13[BMAX][2]; //roots of other functions besides j0 and j1 needed for evaluation of A13 in Appendix.

int static besselsetupxiLPT = 0;
static int xiLPTellmax;
static double besselrootsxiLPT[BMAX][5]; //allow for up to L_8, though I doubt need to keep it.

int fillbesselroots() {
  if(besselsetup == 1) {
    return 0;
    }
  //get spherical bessel function roots.
  double xval;
  double delta = 0.001;
  double fval;
  double fvallast;
  int ello;
  int whichroot;
  double currroot;

  for(ello=0;ello<=1;ello+=1) {
    whichroot=0;
    fvallast=gsl_sf_bessel_jl(ello,0.0);
    xval = delta;
    currroot = 0.0;
    while(currroot < 2500.0) {
      assert(whichroot < BMAX);
      fval=gsl_sf_bessel_jl(ello,xval);
      if((fval > 0.0 && fvallast < 0.0) || (fval < 0.0 && fvallast > 0.0)) {
        besselroots[whichroot][ello] = xval - 0.5*delta;
        if(ello == 0) {
          besselrootsxiLPT[whichroot][ello/2] = xval - 0.5*delta;
          }
        currroot = xval - 0.5*delta;
        whichroot += 1;
        if(whichroot >= BMAX) {
          printf("ahh bessel error! %e\n",currroot);
          exit(1);
          }
        }
      xval += delta;
      fvallast = fval;
      }

    assert(whichroot < BMAX-1);
    }
  besselsetup = 1;
  return 0;
  }

int fillbesselrootsxiLPT(int ellmax) {
  if(besselsetupxiLPT == 1) {
    return 0;
    }
  int ellstart = 0;
  if(besselsetup == 1) {
    ellstart = 2;  //already filled in.
    }

  xiLPTellmax = ellmax;
  assert(xiLPTellmax <=8);
  //get spherical bessel function roots.
  double xval;
  double delta = 0.001;
  double fval;
  double fvallast;
  int ello;
  int whichroot;
  double currroot;

  for(ello=ellstart;ello<=ellmax;ello+=2) {
    whichroot=0;
    fvallast=gsl_sf_bessel_jl(ello,0.0);
    xval = delta;
    currroot = 0.0;
    while(currroot < 2500.0) {
      assert(whichroot < BMAX);
      fval=gsl_sf_bessel_jl(ello,xval);
      if((fval > 0.0 && fvallast < 0.0) || (fval < 0.0 && fvallast > 0.0)) {
        besselrootsxiLPT[whichroot][ello/2] = xval - 0.5*delta;
        currroot = xval - 0.5*delta;
        whichroot += 1;
        if(whichroot >= BMAX) {
          printf("ahh bessel error! %e\n",currroot);
          exit(1);
          }
        }
      xval += delta;
      fvallast = fval;
      }
    assert(whichroot < BMAX-1);
    }
  besselsetupxiLPT = 1;
  return 0;
  }

double getbesselroot(int whichroot, int ell) {
  assert(whichroot <= BMAX-1);
  return(besselroots[whichroot][ell]);
  }

double getbesselrootxiLPT(int whichroot, int ell) {
  assert(whichroot <= BMAX-1);
  assert(ell <= xiLPTellmax);
  return(besselrootsxiLPT[whichroot][ell/2]);
  }

int printbesselrootsxiLPT(int ellmax) {
  int i,j;
  for(i=0;i<BMAX;i++) {
    printf("%d ",i);
    for(j=0;j<=ellmax/2;j++) {
      printf("%e ",getbesselrootxiLPT(i,j*2));
      }
    printf("\n");
    }
  }

double getbesselrootA13(int whichroot, int whichbessel) {
  assert(whichroot < BMAX-1);
  if(whichbessel == 0 || whichbessel == 1) {
    return(besselroots[whichroot][whichbessel]);
    }

  if(whichbessel == 2) {
    return(besselroots[whichroot][1]);
    }
  if(whichbessel == 3 || whichbessel == 4) {
    return(besselrootsA13[whichroot][whichbessel-3]);
    }
  printf("whichbessel = %d out of range in getbesselrootA13\n",whichbessel);
  exit(1);
  }

//additional bessel functions for sigsqr calculation.
//this one is for integration over J(mu^2,kr)
double mybessel(int whichint, double barg) {
  switch(whichint) {
    case(-1):
      return(1.);
    case(0):
      return(gsl_sf_bessel_jl(0,barg));
    case(1):
      if(barg < 0.3) {
        return(1./3.-gsl_pow_2(barg)/30.+gsl_pow_4(barg)/840.); //fractionally accurate to ~1e-7
        }
      else {
        return(gsl_sf_bessel_jl(1,barg)/barg);
        }
    default:
      exit(1);
    }
  }

//this one is for all the random oscillatory functions we need to integrate over in Eqn A13 of the appendix.
double mybesselA13(int whichint, double barg) {
  switch(whichint) {
    case(-1):
      return(1.);
    case(0):
      return(gsl_sf_bessel_jl(0,barg));
    case(1):
      if(barg < 0.3) {
        return(1./3.-gsl_pow_2(barg)/30.+gsl_pow_4(barg)/840.); //fractionally accurate to ~1e-7
        }
      else {
        return(gsl_sf_bessel_jl(1,barg)/barg);
        }
    case(2):
      return(gsl_sf_bessel_jl(1,barg));
    case(3):
      if(barg < 0.3) {
        return(barg*(-1./5.+gsl_pow_2(barg)/42.-gsl_pow_4(barg)/1080.));
        }
      else {
        return(((-6.*barg + gsl_pow_3(barg))*cos(barg) + (6. - 3.*gsl_pow_2(barg))*sin(barg))/gsl_pow_4(barg));
//        return((-6./gsl_pow_3(barg)+1./barg)*cos(barg)+(6./gsl_pow_4(barg)-3./gsl_pow_2(barg)*sin(y));
        }
    case(4):
      if(barg < 0.3) {
        return(barg*(-1./15.+gsl_pow_2(barg)/210.-gsl_pow_4(barg)/7560.));
        }
      else {
        return((3.*barg*cos(barg)+(-3.+gsl_pow_2(barg))*sin(barg))/gsl_pow_4(barg));
        }
    default:
      exit(1);
    }
  }

int fillbesselrootsA13() {
  if(besselsetupA13 == 1) {
    return 0;
    }
  fillbesselroots(); //make sure the other roots are filled as well.
  //get spherical bessel function roots.
  double xval;
  double delta = 0.001;
  double fval;
  double fvallast;
  int whichint;
  int whichroot;
  double currroot;

  for(whichint=3;whichint<=4;whichint+=1) {
    whichroot=0;
    fvallast=mybesselA13(whichint, 0.0); 
    xval = delta;
    currroot = 0.0;
    while(currroot < 2500.0) {
      assert(whichroot < BMAX);
      fval = mybesselA13(whichint, xval); //=gsl_sf_bessel_jl(ello,xval);
      if((fval > 0.0 && fvallast < 0.0) || (fval < 0.0 && fvallast > 0.0)) {
        besselrootsA13[whichroot][whichint-3] = xval - 0.5*delta;
        currroot = xval - 0.5*delta;
        whichroot += 1;
        if(whichroot >= BMAX) {
          printf("ahh bessel error! %e\n",currroot);
          exit(1);
          }
        }
      xval += delta;
      fvallast = fval;
      }
    assert(whichroot < BMAX-1);
    }
  besselsetupA13 = 1;
  return 0;
  }
