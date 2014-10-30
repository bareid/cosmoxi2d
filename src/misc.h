FILE *open_file_write(char *filename);
FILE *open_file_append(char *filename);
FILE *open_file_read(char *filename);
int get_file_length(FILE *ifp, int *headercount, int *linecount);
double **malloc2ddouble(int l1, int l2);
double malloc2dfree(double **m);
double legendrep(int ell, double x);

