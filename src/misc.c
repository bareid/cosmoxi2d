/** BEGIN ROUTINES **/
#include "cosmoxi2d.h"
#include "misc.h"

FILE *open_file_write(char *filename) {
        FILE *ofp;

        if(!(ofp=fopen(filename,"w")))  {
                fprintf(stderr,"Can't open %s\n",filename);
                return NULL;
                }
        else {
                return ofp;
                }
        }

FILE *open_file_append(char *filename) {
        FILE *ofp;

        if(!(ofp=fopen(filename,"a")))  {
                fprintf(stderr,"Can't open %s\n",filename);
                return NULL;
                }
        else {
                return ofp;
                }
        }

FILE *open_file_read(char *filename) {
        FILE *ifp;

        if(!(ifp=fopen(filename,"r")))  {
                fprintf(stderr,"Can't open %s\n",filename);
                return NULL;
                }
        else {
                return ifp;
                }
        }

int get_file_length(FILE *ifp, int *headercount, int *linecount) {
        int i;
        char line[MAXLINELEN];

        *headercount = 0;
        *linecount = 0;

    while(fgets(line,MAXLINELEN,ifp))  {
        if(!feof(ifp))  {
            (*linecount)++;
            }
                else {
                        break;
                        }
                if(strstr(line,"#"))  {
                        (*headercount)++;
                        }
        }
        (void) rewind(ifp);
        }

double **malloc2ddouble(int l1, int l2) {
  double **m;
  int i;
  m = malloc(l1*sizeof(double *));
  m[0] = malloc(l1*l2*sizeof(double));
  for(i=0;i<l1; i++) {
    m[i] = m[0]+i*l2;
    }
  return m;
  }

double malloc2dfree(double **m) {
  if(m) {
    if(m[0]) {
      free(m[0]);
      }
    free(m);
    }
  }

/** END ROUTINES **/

double legendrep(int ell, double x) {

  double xsqr;
  switch(ell) {
    case(0):
  return (1.);
  break;
    case(1):
  return (x);
  break;
    case(2):
  return (-0.5+1.5*x*x);
  break;
    case(3):
  return (-1.5*x+2.5*x*x*x);
  break;
    case(4):
  xsqr = x*x;
  return (0.375 - 3.75*xsqr + 4.375*xsqr*xsqr);
  break;
    case(6):
  xsqr = x*x;
  return (-0.3125+6.5625*xsqr - 19.6875*xsqr*xsqr + 14.4375*xsqr*xsqr*xsqr);
  break;
    default:
  exit(1);
  break;
    }
  }





