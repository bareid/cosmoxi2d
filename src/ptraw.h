#include "cosmoxi2d.h"
#include "linearpk.h"
#include "misc.h"
#include "matsubaraQRandSPT.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int doQRints(double kmin, double kmaxqr1, double kmaxqr2, double dkqr1, double dkqr2);
int doSPTvelints(double kmin, double kmax1, double kmax2, double kmax3, double dk1, double dk2, double dk3);
int getpkrealLPT(double b1L, double b2L, double kminLPT, double kmaxLPT, double dkLPT);
int getpkredLPT(double b1L, double b2L, double fvel, int ellmax, double kminLPT, double kmaxLPT, double dkLPT);
double pkLPT(double k);
double pkredLPT(double k, int ell);
double vinSPTkint(double k);
double sigsqrSPTkint(double k, int whichint);
int getkintvsigrealSPT(double b1L, double b2L, double kmin, double kmax1, double kmax2, double dk1, double dk2);
int readQRints(char *QRfname);
int readSPTvelints(char *SPTfname);
int printQRints(char *QRfname);
//int printQRintsnew(char *QRfname, double kminprint, double kmaxprint, double dkprint);
int printSPTvelints(char *SPTfname);
//int printSPTvelintsnew(char *SPTfname, double kmin, double kmax1, double kmax2, double dk1, double dk2);
int printpkreal(char *pkoutfname);
int printpkredLPT(char *pkoutfname);
int freeptraw();

//temporary checks -- delete access to these eventually.
//double Qkint(double k, int n);
//int printQtilde(double xtest, int n);
//int printRtilde();
