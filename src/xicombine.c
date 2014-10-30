#include "xicombine.h"

static size_t limitint = 20000;
static double epsabs = 1.0e-7;
static double epsrel = 1.0e-7;
static gsl_integration_workspace *workspaceint;

typedef struct {
  int ell;
  double sswitch;
  double sswitchwidth;
  } comboint_params;

double xicombo(double sval, int ellval, double sswitch, double sswitchwidth) {
  if(sval < sswitch - 0.5*sswitchwidth) {
    return(xistreameval(sval,ellval));
    }
  if(sval > sswitch + 0.5*sswitchwidth) {
    return(xistreamLPTeval(sval,ellval));
    }
  double ff = (sval - sswitch + 0.5*sswitchwidth)/sswitchwidth;
  return(xistreameval(sval,ellval)*(1.-ff) + xistreamLPTeval(sval,ellval)*ff);
  }

double xicomboint_integrand(double ss, void *params) {
  comboint_params *p = (comboint_params *) params;
  int ellval = p->ell;
  double sswitch = p->sswitch;
  double sswitchwidth = p->sswitchwidth;
  return(xicombo(ss,ellval,sswitch,sswitchwidth)*ss*ss);
  }

double xicomboint(double binmin, double binmax, int ellval, double sswitch, double sswitchwidth) {
  double comboint_result, comboint_err;
  gsl_function AA;
  int status;
  AA.function = &xicomboint_integrand;
  comboint_params myparams = {ellval,sswitch,sswitchwidth};
  AA.params = &myparams;
  status = gsl_integration_qag(&AA,binmin,binmax, epsabs, epsrel, limitint, GSL_INTEG_GAUSS31, workspaceint, &comboint_result, &comboint_err);

  double vshell = 1./3.*(binmax*binmax*binmax - binmin*binmin*binmin);
  return(comboint_result/vshell);
  }

int printxiscombo(char *outfname, double **sbins, int nsbins, int ellmax,double sswitch, double sswitchwidth,double bias, double fvel, double iso1d, double alphaperp, double alphapar) {

  workspaceint = gsl_integration_workspace_alloc(limitint);

  FILE *ofp;
  int ello,i;
  ofp = open_file_write(outfname);
  fprintf(ofp,"b = %f\n",bias);
  fprintf(ofp,"f = %f\n",fvel);
  fprintf(ofp,"iso1d = %f\n",iso1d);
  fprintf(ofp,"alphaperp = %f\n",alphaperp);
  fprintf(ofp,"alphapar = %f\n",alphapar);
  for(i=0;i<nsbins;i++) {
    fprintf(ofp,"%e ",sbins[i][0]);
    for(ello=0;ello<=ellmax;ello+=2) {
      fprintf(ofp,"%e ",xicomboint(sbins[i][1],sbins[i][2],ello,sswitch,sswitchwidth));
      }
    fprintf(ofp,"\n");
    }
  gsl_integration_workspace_free(workspaceint);
  return 0;
  }


int printxiscombonoint(char *outfname, double **sbins, int nsbins, int ellmax,double sswitch, double sswitchwidth,double bias, double fvel, double iso1d, double alphaperp, double alphapar) {
  printf("hey beth inside\n");
  FILE *ofp;
  int ello,i;
  ofp = open_file_write(outfname);
  fprintf(ofp,"b = %f\n",bias);
  fprintf(ofp,"f = %f\n",fvel);
  fprintf(ofp,"iso1d = %f\n",iso1d);
  fprintf(ofp,"alphaperp = %f\n",alphaperp);
  fprintf(ofp,"alphapar = %f\n",alphapar);
  for(i=0;i<nsbins;i++) {
    fprintf(ofp,"%e ",sbins[i][0]);
    for(ello=0;ello<=ellmax;ello+=2) {
      fprintf(ofp,"%e ",xicombo(sbins[i][0],ello,sswitch,sswitchwidth));
      }
    fprintf(ofp,"\n");
    }
  return 0;
  }


