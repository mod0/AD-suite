#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

extern double gam, gm1, cfl, eps, mach, alpha;
void input(int maxnode, int maxcell, int maxedge, int& nnode, int& ncell,
	   int& nedge, double **x, double** q, int** cell, int** edge, 
	   int** ecell, int* boun);

int time_cell(double* x1, double* x2, double* x3, double* x4, double* q,
	      double &adt);

int flux_face(double* x1, double* x2, double* q1, double* q2, double adt1, 
	      double adt2, double* res1, double* res2);

int flux_wall(double* x1, double* x2, double* q, double* res);

int lift_wall(double* x1, double* x2, double* q, double &lift);

double max(double arg1,double arg2);

double min(double arg1, double arg2);

double abs(double arg1);

void output(int ncell, double** q);
#endif
