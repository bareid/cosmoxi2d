#include "cosmoxi2d.h"
#include "misc.h"
#include "linearpk.h"
#include "ptraw.h"
#include "peakbackground.h"
#include "pktoxi.h"
#include "pktoxilin.h"
#include "pktoxiLPT.h"
#include "pktovin.h"
#include "pktosig.h"
#include "xistream.h"
#include "xistreamLPT.h"
#include "xicombine.h"

//struct for holding input parameters.
typedef struct {
  int ptopt;
  char pklinfname[MAXLINELEN];
  double pkhval;
  double pkampscale;
  char qrintsfname[MAXLINELEN];
  char sptvelfname[MAXLINELEN];
  char xi024outfname[MAXLINELEN];
  double bias;  //galaxy bias on large scales.
  double iso1d; //additional incoherent isotropic 1d velocity dispersion.
  double fvel;  //amplitude of quasilinear velocity field (std "f" parameter = d ln D/d ln a)
  double alphaperp; //stretch in the perpendicular direction
  double alphapar; //stretch in the LOS direction.
  double deltasfine; //used internally to avoid interpolation errors; linear spacing in s
  char databinfile[MAXLINELEN];
  int ellmax; //how many Legendre moments do you need
  } runparams;

const int nrunparams = 15;

//these are filled in in getbins()
double **sbins;  //2d array [nsbins][3], each rbin contains center and two edges of the bin for clarity.
int nsbins;
double hfid;  //not used now, may want later.

int getbins(char *databinfile) {
  FILE *ifp;
  char line[MAXLINELEN];
  char *r;
  int i;
  int headercount, linecount;
  
  ifp = open_file_read(databinfile);
  get_file_length(ifp,&headercount,&linecount);
  assert(headercount == 1); //corresponds to hfid
  fgets(line,MAXLINELEN,ifp);
  assert(strstr(line,"#"));
  r = strchr(line,'=');
  assert(r);
  hfid = ((double) atof(r+1));
  nsbins = linecount - headercount;
  sbins = malloc2ddouble(nsbins,3);
  for(i=0;i<nsbins;i++) {
    fscanf(ifp,"%le %le %le\n",&(sbins[i][0]),&(sbins[i][1]),&(sbins[i][2]));
    }
  fclose(ifp);
  }

runparams parserunfile(char *infilename) {
  runparams myparams;
  FILE *ifpparams;
  ifpparams = open_file_read(infilename);
  char line[MAXLINELEN];
  char *r;
  int mycount = 0;

  int gotlist[nrunparams];
  int i;

  for(i=0;i<nrunparams;i++) {
    gotlist[i] = 0;
    }

  while(fgets(line,MAXLINELEN,ifpparams))  {
    if(strstr(line,"#"))  {
      continue;
      }
    if(!feof(ifpparams))  {
      if(strstr(line,"ptcalcopt")) {
        r = strchr(line,'=');
        assert(r);
        myparams.ptopt = atoi(r+1);
        assert(myparams.ptopt >= 0 && myparams.ptopt <=3);
        gotlist[0] = 1;
        continue;
        }
      if(strstr(line,"pklinfname")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.pklinfname);
        gotlist[1] = 1;
        continue;
        }
      if(strstr(line,"pkhval")) {
        r = strchr(line,'=');
        assert(r);
        myparams.pkhval = atof(r+1);
        gotlist[2] = 1;
        continue;
        }
      if(strstr(line,"pkampscale")) {
        r = strchr(line,'=');
        assert(r);
        myparams.pkampscale = atof(r+1);
        gotlist[14] = 1;
        continue;
        }
      if(strstr(line,"qrintsfname")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.qrintsfname);
        gotlist[3] = 1;
        continue;
        }
      if(strstr(line,"sptvelfname")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.sptvelfname);
        gotlist[4] = 1;
        continue;
        }
      if(strstr(line,"xi024outfname")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.xi024outfname);
        gotlist[5] = 1;
        continue;
        }
      if(strstr(line,"bias")) {
        r = strchr(line,'=');
        assert(r);
        myparams.bias = atof(r+1);
        gotlist[6] = 1;
        continue;
        }
      if(strstr(line,"isotropicdisp1d")) {
        r = strchr(line,'=');
        assert(r);
        myparams.iso1d = atof(r+1);
        gotlist[7] = 1;
        continue;
        }
      if(strstr(line,"fvel")) {
        r = strchr(line,'=');
        assert(r);
        myparams.fvel = atof(r+1);
        gotlist[8] = 1;
        continue;
        }
      if(strstr(line,"alphaperp")) {
        r = strchr(line,'=');
        assert(r);
        myparams.alphaperp = atof(r+1);
        gotlist[9] = 1;
        continue;
        }
      if(strstr(line,"alphapar")) {
        r = strchr(line,'=');
        assert(r);
        myparams.alphapar = atof(r+1);
        gotlist[10] = 1;
        continue;
        }
      if(strstr(line,"deltasfine")) {
        r = strchr(line,'=');
        assert(r);
        myparams.deltasfine = atof(r+1);
        gotlist[11] = 1;
        continue;
        }
      if(strstr(line,"databinfile")) {
        r = strchr(line,'=');
        assert(r);
        sscanf(r+1,"%s\n",myparams.databinfile);
        gotlist[12] = 1;
        continue;
        }
      if(strstr(line,"ellmax")) {
        r = strchr(line,'=');
        assert(r);
        myparams.ellmax = atoi(r+1);
        gotlist[13] = 1;
        continue;
        }
      }
    else {
      break;
      }
    }
  (void) rewind(ifpparams);
  fclose(ifpparams);

  // do some checks on the input parameters.
  for(i=0;i<nrunparams;i++) {
    //don't need file names for qrints and sptvel if you're going to compute them.
    if(i==3 && gotlist[i] == 0) {
      assert(myparams.ptopt == 2 || myparams.ptopt == 3);
      sprintf(myparams.qrintsfname,"qrintsnew");
      continue;
      }
    if(i==4 && gotlist[i] == 0) {
      assert(myparams.ptopt == 1 || myparams.ptopt == 3);
      sprintf(myparams.sptvelfname,"sptvelintsnew");
      continue;
      }
    if(gotlist[i] == 0) {
      printf("missing %d entry of params file\n",i+1);
      }
    assert(gotlist[i] == 1);
    }
  return myparams;
  }

int main(int argc, char *argv[]) {
  if(argc != 2) {
    fprintf(stderr,"Usage: ./cosmoxi2d paramfile\n");
    exit(1);
    }

  time_t time1, time2;
  (void) time(&time1);

  char paramfile[MAXLINELEN];
  sprintf(paramfile,"%s",argv[1]);

  runparams myparams = parserunfile(paramfile);
  getbins(myparams.databinfile);
  double bias = myparams.bias;
  int ptopt = myparams.ptopt;
  double b1L, b2L;

  getpkspline(myparams.pklinfname, myparams.pkhval, myparams.pkampscale); //set up linear pk spline; internal units are Mpc^-1.  Assuming input P(k) from camb is in units h Mpc^-1, we convert to Mpc^-1.  If the input P(k) is in Mpc^-1, then pass 1.0 for pkhval.
  //integration variables we hold fixed, but probably want to adjust.
  //these were tweaked in h Mpc^-1 with h=0.7, so adjust accordingly.
  double kminpkLPT = 0.001*0.7;
  double kmaxpkLPT = 1.0*0.7;
  double dkLPT = 0.001*0.7;
  double dk2vel = 0.05*0.7;

  double kmin = 0.001*0.7;
  double kmaxqr1 = 0.3*0.7;
  double kmaxqr2 = 1.0*0.7;
  double dkqr1 = 0.005*0.7*3.;
  double dkqr2 = 0.03*0.7*3.;

  double kmaxv1 = 0.4*0.7;
  double kmaxv2 = 1.0*0.7;
  double kmaxv3 = 10.0*0.7;
  double dkv1 = 0.005*0.7*3.;
  double dkv2 = 0.03*0.7*3.;
  double dkv3 = 0.25*0.7*3.;

  double kmaxxi= kmaxv3; //THIS IS USED TO INTEGRATE XILIN.

  double sswitch = 70./0.7;
  double sswitchwidth = 5./0.7;

  //extract these from redge/rcentre files.
  double smin = sbins[0][1];
  double smax = sbins[nsbins-1][2];


  //go beyond where we want xi024 so we know the integrand of the streaming model well.
  double rmin = smin - floor(smin/myparams.deltasfine)*myparams.deltasfine;
  assert(rmin > 0.);
  double rmax = sswitch+0.5*sswitchwidth + 60./0.7; //really generous overlap.
  double dr = myparams.deltasfine;  //can adjust this later.

  double rminLPT = max(sswitch-0.5*sswitchwidth - sqrt(myparams.iso1d)*10.,rmin);
  //original.
  //double rmaxLPT = smax + sqrt(myparams.iso1d)*10.;
  double rmaxLPT = smax + sqrt(myparams.iso1d)*12.;
  double drLPT = dr; //if i want these unequal, need to rewrite xicombo.

  if(ptopt == 1 || ptopt == 3) {
    doSPTvelints(kmin,kmaxv1,kmaxv2,kmaxv3,dkv1,dkv2,dkv3);
    printSPTvelints(myparams.sptvelfname);
    }
  else {
    readSPTvelints(myparams.sptvelfname);
    }

  if(ptopt == 2 || ptopt == 3) {
    doQRints(kmin,kmaxqr1,kmaxqr2,dkqr1,dkqr2);
    (void) time(&time2);
    #ifdef VERBOSE
    printf("doQRints done.  time since start = %e\n",difftime(time2,time1));
    #endif
    printQRints(myparams.qrintsfname);
    }
  else {
    readQRints(myparams.qrintsfname);
    }
  getbL(bias,&b1L,&b2L);
  #ifdef VERBOSE
  printf("bias, b1L, b2L: %e %e %e\n",bias,b1L,b2L);
  #endif

  getpkrealLPT(b1L,b2L,kminpkLPT,kmaxpkLPT,dkLPT);
  getpkredLPT(b1L,b2L,myparams.fvel,myparams.ellmax,kminpkLPT,kmaxpkLPT,dkLPT);

  /*
  char pkLPTredout[MAXLINELEN];
  sprintf(pkLPTredout,"pkLPTredb%.3f.dat",bias);
  printpkredLPT(pkLPTredout);
  */

  (void) time(&time2);
  #ifdef VERBOSE
  printf("pkredLPT done.  time since start = %e\n",difftime(time2,time1));
  #endif
  getxiLPTmoments(kmaxpkLPT,rminLPT,rmaxLPT,drLPT,myparams.ellmax);
  /*
  char xiLPTredoutfname[MAXLINELEN];
  sprintf(xiLPTredoutfname,"xiLPTredb%.3f.dat",bias);
  printxiLPTmoments(xiLPTredoutfname);
  */
  (void) time(&time2);
  #ifdef VERBOSE
  printf("xiLPTmoments done.  time since start = %e\n",difftime(time2,time1));
  #endif
  dostreamxiLPT(sswitch-0.5*sswitchwidth,smax,myparams.deltasfine,myparams.iso1d,myparams.ellmax, myparams.alphaperp, myparams.alphapar);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("streamLPT done.  time since start = %e\n",difftime(time2,time1));
  #endif
  /*
  char xiLPTstreamoutfname[MAXLINELEN];
  sprintf(xiLPTstreamoutfname,"xiLPTstreamcompareredb%.3f.dat",bias);
  printxisLPT(xiLPTstreamoutfname,myparams.bias,myparams.fvel,myparams.iso1d,myparams.alphaperp, myparams.alphapar);
  */
  getkintvsigrealSPT(b1L,b2L,kminpkLPT,kmaxpkLPT,kmaxv3,dkLPT,dk2vel);
  //printSPTvelintsnew(tmpfname,kminpkLPT,kmaxpkLPT,kmaxv3,dkLPT,dk2vel);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getkintvsigrealSPT done.  time since start = %e\n",difftime(time2,time1));
  #endif
  getxireal(kmaxpkLPT,rmin,rmax,dr);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getxireal done.  time since start = %e\n",difftime(time2,time1));
  #endif
  getxilini(kmaxxi,rmin,rmax,dr,bias);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getxilini done.  time since start = %e\n",difftime(time2,time1));
  #endif
  //printxireal(xioutfname);
  //printxilin(xilinoutfname);
  getvinreal(kmaxv3,rmin,rmax,dr);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getvinreal done.  time since start = %e\n",difftime(time2,time1));
  #endif
  getvinlinreal(kmaxv3,rmin,rmax,dr,b1L);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getvinlinreal done.  time since start = %e\n",difftime(time2,time1));
  #endif
  //printvinreal(vinoutfname);
  //printvinlinreal(vinlinoutfname);
  getsigsqrreal(kmaxv3,rmin,rmax,dr,b1L);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("getsigsqrreal done.  time since start = %e\n",difftime(time2,time1));
  #endif
  //printsigsqr(sigsqroutfname);
  (void) time(&time2);
  #ifdef VERBOSE
  printf("real space xi,vin,sigsqr done.  time since start = %e\n",difftime(time2,time1));
  #endif

  //do streaming model convolution.
  dostream(smin,sswitch+0.5*sswitchwidth,myparams.deltasfine,myparams.iso1d,myparams.fvel,myparams.ellmax, myparams.alphaperp, myparams.alphapar); 
  (void) time(&time2);
  #ifdef VERBOSE
  printf("stream done.  time since start = %e\n",difftime(time2,time1));
  #endif
  //want to spit out all the relevant parameters used to compute this xi024.
  /*
  char xistreamout[MAXLINELEN];
  sprintf(xistreamout,"xistreamoutb%.3f.dat",bias);
  printxis(xistreamout,myparams.bias,myparams.fvel,myparams.iso1d,myparams.alphaperp, myparams.alphapar);
  */

  //testing.
  /*
  char xioutwint[MAXLINELEN];
  sprintf(xioutwint,"xitestoutwintb%.3f.dat",bias);
  printxiscombo(xioutwint, myparams.smin, myparams.smax,myparams.deltasdata, myparams.ellmax, sswitch,sswitchwidth,myparams.bias,myparams.fvel,myparams.iso1d,myparams.alphaperp, myparams.alphapar);

  char xioutnoint[MAXLINELEN];
  sprintf(xioutnoint,"xitestoutnointb%.3f.dat",bias);
  printxiscombonoint(xioutnoint, myparams.smin, myparams.smax,myparams.deltasdata, myparams.ellmax, sswitch,sswitchwidth,myparams.bias,myparams.fvel,myparams.iso1d,myparams.alphaperp, myparams.alphapar);
  //end testing.
  */
  printxiscombo(myparams.xi024outfname, sbins, nsbins, myparams.ellmax, sswitch,sswitchwidth,myparams.bias,myparams.fvel,myparams.iso1d,myparams.alphaperp, myparams.alphapar);
  //cleanup
  pkfree();
  freeptraw();
  pktoxifree();
  pktoxiLPTfree();
  pktovinfree();
  pktosigsqrfree();
  pktoxilinfree();
  streamfree();
  streamLPTfree();
  malloc2dfree(sbins);
  return 0;
  }
