#include "peakbackground.h"

/******
This code is used to solve for b1L and b2L in terms of input Eulerian bias using the Sheth-Tormen mass function.
These expressions are from v4 on the arXiv:0807.1733v4; note there was an erratum that updated Eqns 57-58 (PRD 78 109901), which give <F'> and <F''>
The algorithm is to solve for the qnu2 that gives <F'> = bEulerian - 1, and plug in to get <F''>
*******/

//some constants associated with ST mass fxn.
const double delc = 1.686;
const double p = 0.3;

struct solve_params {
  double bfit;
  };

double Fp(double qnu2) {
  return((qnu2 - 1. + 2.*p/(1.+pow(qnu2,p)))/delc);
  }

double Fpp(double qnu2) {
  return((qnu2*qnu2 - 3.*qnu2 + 2.*p*(2.*qnu2+2.*p - 1.)/(1.+pow(qnu2,p)))/(delc*delc));
  }

double fsolve(double x, void *params) {
  struct solve_params *p = (struct solve_params *) params;
  double bfit = p->bfit;
  return(Fp(x)-bfit);
  }
double fsolve_deriv(double x, void *params) {
  return((1. - 2.*p*p*pow(x,p-1.)/gsl_pow_2(1.+pow(x,p)))/delc);
  }

void fsolve_fdf(double x, void *params, double *y, double *dy) {
  *y = fsolve(x,params);
  *dy = fsolve_deriv(x,params);
  }

//copying from http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html
int getbL(double bias, double *b1L, double *b2L) {

  if(bias-1. < (-1.+2.*p)/delc) {
    printf("bias smaller than minimum: %e %e\n",bias,(-1.+2.*p)/delc);
    exit(1);
    }

  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x_lo = 0.0, x_hi = 1.1*delc*bias+1.;
  double x = x_hi;
  double x0;
  assert(Fp(x_lo) < bias-1. && Fp(x_hi) > bias-1.);
  gsl_function_fdf FDF;

  FDF.f = &fsolve;
  FDF.df = &fsolve_deriv;
  FDF.fdf = &fsolve_fdf;
  struct solve_params myp = {bias-1.};
  FDF.params = &myp;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x_hi);

  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1.0e-6);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);

  (*b1L) = Fp(x);
  if(!(fabs((*b1L) + 1. - bias) < 2.0e-6)) {
    printf("bias-1 + b1L: %e %e %e\n",bias-1+(*b1L), bias,(*b1L));
    }
  assert(fabs((*b1L) + 1. - bias) < 2.0e-6);
  (*b2L) = Fpp(x);
  assert(status == 0);
  return status;
  }


