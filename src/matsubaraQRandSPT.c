#include "matsubaraQRandSPT.h"

/******
Qtilde(1-13) corresponds to Qtilde in Matsubara 2008b (PRD 78, 083519), Eqn A39-A46.
Rtilde(1-2) corresponds to Rtilde in Matsubara 2008b (PRD 78, 083519), Eqn A47-A48.

We've tacked on relevant integrals from SPT Pdt and Ptt that are needed for the velocity calculations as well.

Q(14) corresponds to P^(2,2)_dt (eqn A26 of Carlson, White, Padmanabhan 2008 (PRD 80 043531)
Q(15) corresponds to P^(2,2)_dt (eqn A24 of Carlson, White, Padmanabhan 2008 (PRD 80 043531)

*****/

double Qtilde(int n, double r, double x) {
  double rsqr, xsqr,denom;
  double dx, dr;
  switch(n) {
    case(1):
  if(r == 1. && x == 1.) {
    return(1.);
    }
  if(x == 1.) {
    return(0.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  if(denom < 0.1) {
    //rescale in terms of dr = 1-r, dx = 1-x; dx^2 in numerator and denominator cancel.
    dx = 1.-x;
    dr = 1.-r;
    return(rsqr*gsl_pow_2(-1.+dx*0.5)/gsl_pow_2(1.-dr+0.5*dr*dr/dx));
    }
  return(rsqr*gsl_pow_2(1.-x*x)/gsl_pow_2(denom));
     case(2):
  if(r == 1. && x == 1.) { //no well-defined limit.
    return(0.);
    }
  if(x == 1.) {
    return(0.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  if(denom < 0.1) {
    //rescale in terms of dr = 1-r, dx = 1-x; dx^2 in numerator and denominator cancel.
    dx = 1.-x;
    dr = 1.-r;
    return(0.5*(1.-dx*0.5)*(1.-dr-dx+dr*dx)*(1.+dr/dx-dr)/gsl_pow_2(1.-dr+0.5*dr*dr/dx));
    }
  return((1.-x*x)*r*x*(1.-r*x)/gsl_pow_2(denom));
    case(3):
  //Q3 diverges as 1/(1-r)^2 near r=1 at x=1.  We'll do r integrals and then put x=1 as a singular point.
  if(r == 1. && x == 1.) {
    return(0.25);  //this is true at fixed r=1, limit as x goes to 1.
    }
  if(x == 1.) {
    return(1./gsl_pow_2(1.-r));
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(x*x*gsl_pow_2(1.-r*x)/gsl_pow_2(denom));
    case(4):
    //diverges as 1/(1-x) for r=1 near x=1.
  if(x == 1.) {
    return 0.;
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return((1.-x*x)/gsl_pow_2(denom));
    case(5):
  if(r == 1. && x == 1.) {
    return(1.);
    }
  if(x == 1.) {
    return(0.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(r*x*(1.-x*x)/denom);

    case(6):
  if(r == 1. && x == 1.) {
    return(-2.);
    }
  if(x == 1.) {
    return(0.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return((1.-3.*r*x)*(1.-x*x)/denom);

    case(7):
  if(r == 1. && x == 1.) { //at x=1, as r -> 1, diverges as 1/(1-r)
    return(0.5);
    }
  if(x == 1.) {
    return(1./(1.-r));
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(x*x*(1.-r*x)/denom);
  exit(1);

    case(8):
  if(r == 1. && x == 1.) {
    return(1.);
    }
  if(x == 1.) {
    return(0.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(rsqr*(1.-x*x)/denom);

    case(9):
  if(r == 1. && x == 1.) { //at x=1, as r -> 1, diverges as 1/(1-r)
    return(0.5);
    }
  if(x == 1.) {
    return(r/(1.-r));
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(r*x*(1.-r*x)/denom);

    case(10):
  return(1.-x*x);
    case(11):
  return(x*x);
    case(12):
  return(r*x);
    case(13):
  return(r*r);
    case(14): //A26 of Carlson, Pdt^(2,2)
  if(r == 1. && x == 1.) { //avoid singularity.
    return(0.);
    }
  if(x == 1.) {  //at x=1, diverges!!
    return(49./gsl_pow_2(1.-r)/98.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return((3.*r+7.*x-10.*r*x*x)*(7.*x - r - 6.*r*x*x)/gsl_pow_2(denom)/98.);

    case(15): //A24 of Carlson, Ptt^(2,2)
  if(r == 1. && x == 1.) {
    return(0.);
    }
  if(x == 1.) {  //at x=1, diverges!!
    return(49./gsl_pow_2(1.-r)/98.);
    }
  rsqr = gsl_pow_2(r);
  denom = (1.+rsqr-2.*r*x);
  return(gsl_pow_2(7.*x - r - 6.*r*x*x)/gsl_pow_2(denom)/98.);

  //velocity integrals.
    case(16):  //Eqn A6 of paper for infall velocity.
  if(r == 1. && x == 1.) {
    return(0.);
    }
  if(x == 1.) {  //at x=1, diverges!!
    return(0.5/(1.-r)); 
    }
  rsqr = gsl_pow_2(r);
  xsqr = gsl_pow_2(x);
  denom = (1.+rsqr-2.*r*x);
  return((3.*r*x-10.*r*x*xsqr+7.*xsqr)/(14.*denom));
    case(17):  //Eqn A14 of paper for sigma2 Bdtt term.
  if(r == 1. && x == 1.) {
    return(0.);
    }
  if(x == 1.) {  //at x=1, diverges!!
    return(7./(r-1.)); 
    }
  rsqr = gsl_pow_2(r);
  xsqr = gsl_pow_2(x);
  denom = (1.+rsqr-2.*r*x);
  return((r*x + 6.*r*x*xsqr - 7.*xsqr)/(14.*denom));
     default:
  exit(1);    
  } //end switch
  }

//copied from SPTinfallvel.c
//y function in Eqn A5 of Reid and White appendix.
double myyfuncvin(double y) {
  double ysqr = gsl_pow_2(y);
  double ysqrinv;
  double atanhdiff;
  double myresult;
  if(y < 0.5) {
    myresult = -1./3.+4./7.*ysqr-12./35.*gsl_pow_2(ysqr)+1./420.*(27.*gsl_pow_3(ysqr)-12.*gsl_pow_4(ysqr)+9.*gsl_pow_5(ysqr));
  if(y > 0.001) { //to avoid problems with 1/y
    atanhdiff = atanh(2.*y/(1.+ysqr)) - 2.*y*(1.+ysqr/3.+gsl_pow_2(ysqr)/5.);
    myresult += 9./168.*gsl_pow_3(-1.+ysqr)*(atanhdiff)/y;
    }
    return myresult;
    }
  else {
    if(y > 2.0) {
      ysqrinv = 1./ysqr;
      myresult = 1./105. - 12./245.*ysqrinv -17./980.*gsl_pow_2(ysqrinv) + 6./245.*gsl_pow_3(ysqrinv) - 3./196.*gsl_pow_4(ysqrinv);
      atanhdiff = atanh(2.*y/(1.+ysqr)) - 2./y*(1.+ysqrinv/3.+gsl_pow_2(ysqrinv)/5.+gsl_pow_3(ysqrinv)/7.);
      if(y<30.) { //after 30, it becomes a 1.5e-6 correction; this is causing us a convergence problem anyway.
        myresult += 9./168.*gsl_pow_3(-1.+ysqr)*(atanhdiff)/y;
        }
      return myresult;
      }
    else {
      //just evaluate the function.
      myresult = 2.*y*(9.+52.*ysqr-9.*gsl_pow_2(ysqr)) - 56.*y*(1.+ysqr);
      if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
        myresult += 9.*gsl_pow_3(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
        }
      myresult = myresult/168./y;
      return myresult;
      }
    }
  }

double myyfuncsigsqr(double y, int whichyint) {
  double ysqr = gsl_pow_2(y);
  double ysqrinv;
  double atanhdiff;
  double myresult;
  switch(whichyint) {
    case(6):
  //this time, just did real expansion instaed of arctanhdiff
  if(y < 0.5) {
    myresult = 8./35.*ysqr-24./245.*gsl_pow_2(ysqr)+8./735.*(gsl_pow_3(ysqr)+gsl_pow_4(ysqr)/11.+3./143.*gsl_pow_5(ysqr));
    return myresult;
    }
  else {
    if(y > 2.0) {
      ysqrinv = 1./ysqr;
      myresult = 8./35.-24./245.*ysqrinv+8./735.*(gsl_pow_2(ysqrinv)+gsl_pow_3(ysqrinv)/11.+3./143.*gsl_pow_4(ysqrinv)+1./143.*gsl_pow_5(ysqrinv));
      return myresult;
      }
    else {
      //just evaluate the function.
      myresult = 2.*y*(-9.+33.*ysqr+33.*gsl_pow_2(ysqr)-9.*gsl_pow_3(ysqr));
      if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
        myresult += 9.*gsl_pow_4(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
        }
      myresult = myresult/(672.*ysqr*y);
      return(myresult);
      }
    }

//SWITCH SIGN ON THIS ONE FOR CONVENIENCE!  IT STAYS POSITIVE.
    case(7):
  if(y < 0.5) {
    myresult = -1./3.+12./35.*ysqr-12./49.*gsl_pow_2(ysqr)+4./105.*gsl_pow_3(ysqr)+12./2695.*gsl_pow_4(ysqr)+4./3185.*gsl_pow_5(ysqr);  //agrees with true result up to 2e-7
    return(-myresult);
    }
  else {
    if(y > 2.0) {
      ysqrinv = 1./ysqr;
      myresult = -23./105. + 12./245.*ysqrinv -4./245.*gsl_pow_2(ysqrinv) - 4./1617.*gsl_pow_3(ysqrinv) - 4./5005.*gsl_pow_4(ysqrinv)-12./35035.*gsl_pow_5(ysqrinv);  //agrees at 5e-8 level
      return(-myresult);
      }
    else {
      //just evaluate the function.
      myresult = 2.*y*(9.-109.*ysqr+63.*gsl_pow_2(ysqr)-27.*gsl_pow_3(ysqr));
      if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.  leaving this out is ~1e-5 difference at 0.001
        myresult += 9.*gsl_pow_3(-1.+ysqr)*(1.+3.*ysqr)*atanh(2.*y/(1.+ysqr));
        }
      myresult = myresult/(672.*ysqr*y);
      return(-myresult);
      }
    }
    case(8): //this is the sume of d1v2v1xpr + d2v1v1xpr, so that the integrals are finite.
  if(y < 0.5) {
    myresult = -1./3. + 4./7.*ysqr-12./35.*gsl_pow_2(ysqr)+12./245.*gsl_pow_3(ysqr)+4./735.*gsl_pow_4(ysqr)+4./2695.*gsl_pow_5(ysqr);
    return(-myresult);
    }
  else {
    if(y > 2.0) {
      ysqrinv = 1./ysqr;
      myresult = 1./105.-12./245.*ysqrinv-4./735.*gsl_pow_2(ysqrinv)-4./2695.*gsl_pow_3(ysqrinv)-4./7007.*gsl_pow_4(ysqrinv)-4./15015.*gsl_pow_5(ysqrinv);
      return(-myresult);
      }
    else {
      //just evaluate the function.
      myresult = -2.*y*(19.-24.*ysqr+9.*gsl_pow_2(ysqr));
      if(fabs(y-1) > 0.001) { //stay away from singularity at y=1; evaluates to 0 anyway for this term.
        myresult += 9.*gsl_pow_3(-1.+ysqr)*atanh(2.*y/(1.+ysqr));
        }
      myresult = myresult/(168.*y);
      return(-myresult);
      }
    }
    default:
  exit(1);
    }
  }

double Rtilde(int n, double r) {
  double rsqr = r*r;

  switch(n) {
    case(1):
  if(r < 0.05) { //fractional accuracy is 3e-7
    return(16.*(rsqr/15.-rsqr*rsqr/35.));
    }
  if(r > 0.985 && r < 1.015) { //fractional accuracy is 3e-7
    return(2./3.*(1.+(r-1)-gsl_pow_2(r-1)));
    }
  if(r > 5.) { //fractional accuracy is 3e-7
    return(16.*(1./15.-1./35./rsqr+1./315./gsl_pow_2(rsqr)));
    }
  return(-(1.+rsqr)/(24.*rsqr)*(3.-14.*rsqr+3.*rsqr*rsqr) + gsl_pow_4(rsqr-1.)/(16.*r*rsqr)*log(fabs((1.+r)/(1.-r))));
    case(2):
  if(r < 0.045) { //fractional accuracy 1.0e-6
    return(4./15.*rsqr - 12./35.*rsqr*rsqr);
    }
  if(r > 0.994 && r < 1.006) { //fractional accuracy 0.0005, absolute accuracy 1.0e-6
    return((1.-r)/3.-gsl_pow_2(r-1.)/6.);
    }
  if(r > 7.) { //fractional accuracy is 3e-7
    return(-4./15.+12./35./rsqr-4./63./gsl_pow_2(rsqr));
    }
  return((1.-rsqr)/(24.*rsqr)*(3.-2.*rsqr+3.*rsqr*rsqr) + (1.+rsqr)*gsl_pow_3(rsqr-1.)/(16.*r*rsqr)*log(fabs((1.+r)/(1.-r))));

    case(3): //Pdt^(1,3) Eqn A25
  if(r < 0.1) {  //fractional accuracy 1e-7
    return((-168.+416./5.*rsqr-2976./35.*rsqr*rsqr)/252.);
    }
  if(fabs(r-1.) < 0.005) { //fractionally accurate 1e-6
    return((-152.-56.*(r-1.)-52.*gsl_pow_2(r-1.))/252.);
    }
  if(r > 10.) { //fractionally accurate 1.0e-8
    return((-200.+2208./35./rsqr-1312./105./gsl_pow_2(rsqr))/252.);
    }
  return((24./rsqr-202.+56.*rsqr-30.*rsqr*rsqr+3./(rsqr*r)*gsl_pow_3(rsqr-1.)*(5.*rsqr+4.)*log(fabs((1.+r)/(1.-r))))/252.);

    case(4): //Ptt^(1,3) Eqn A23
  if(r < 0.1) {  //fractional accuracy 1e-7
    return((-56.-32./5.*rsqr-96./7.*rsqr*rsqr)/84.);
    }
  if(fabs(r-1.) < 0.005) { //fractionally accurate 1e-6
    return((-72.-40.*(r-1.)+4.*gsl_pow_2(r-1.))/84.);
    }
  if(r > 10.) { //fractionally accurate 1.0e-8
    return((-504./5.+1248./35./rsqr-608./105./gsl_pow_2(rsqr))/84.);
    }
  return((12./rsqr-82.+4.*rsqr-6.*rsqr*rsqr+3./(rsqr*r)*gsl_pow_3(rsqr-1.)*(rsqr+2.)*log(fabs((1.+r)/(1.-r))))/84.);

    case(5): //Eqn A5
  return(myyfuncvin(r));
    case(6):  //Eqn A12 top one.
  return(myyfuncsigsqr(r,6));
    case(7): //Eqn A12 bottom one.
  return(myyfuncsigsqr(r,7));
    case(8): //Eqn A15
  return(myyfuncsigsqr(r,8));
    } //end switch. 
  exit(1);
  }
