#include <iostream>
using namespace std;

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "const.h"


double arg1,arg2;



/********************************************************************************************************/


int time_cell(double* x1, double* x2, double* x3, double* x4, double* q,
	      double& adt)
	      
{	      

//-----------------------compute area/timestep for an individual cell-----------------------------------------	      

  int i; 
  double ri, u, v, p, c;
  double *dx = new double[4];
  double *dy = new double[4];


	dx[0] = x2[0] - x1[0];
	dx[1] = x3[0] - x2[0];
	dx[2] = x4[0] - x3[0];
	dx[3] = x1[0] - x4[0];

	dy[0] = x2[1] - x2[1];
	dy[1] = x3[1] - x2[1];
	dy[2] = x4[1] - x3[1];
	dy[3] = x1[1] - x4[1];
	
	ri = 1.0/q[0];
	u = ri*q[1];
	v =  ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));
	c = sqrt(gam*ri*p);
	
	adt = 0;

	for (i = 0; i < 4; i++) 
	{
		adt = adt + fabs(u*dy[i] - v*dx[i]) + c*sqrt(dx[i]*dx[i] 
		      + dy[i]*dy[i]); 
	} 

	adt = adt/cfl;

delete[] dx;
delete[] dy;

return 0;

} //end


/*********************************************************************************************************/

int flux_face(double* x1, double* x2, double* q1, double* q2, double adt1, 
	      double adt2, double* res1, double* res2)
	      
{	      
	      
//------------------------compute flux through an individual face-----------------------------------------	      


int n, r; 
double dx, dy, mu, r1i, u1, v1, p1, vol1, r2i, u2, v2, p2, vol2;
double *f = new double[4];

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];
	mu = 0.5*(adt1 + adt2)*eps;

	r1i = 1.0/q1[0];
	u1 = r1i*q1[1];
	v1 = r1i*q1[2];
	p1 = gm1*(q1[3] - 0.5*r1i*(q1[1]*q1[1] + q1[2]*q1[2]));
	vol1 = u1*dy - v1*dx;
	
	r2i = 1.0/q2[0];
	u2 = r2i*q2[1];
	v2 = r2i*q2[2];
	p2 = gm1*(q2[3] - 0.5*r2i*(q2[2]*q2[2] + q2[3]*q2[3]));
	vol2 = u2*dy - v2*dx;
	
	f[0] = 0.5*(vol1*q1[0]         + vol2*q2[0]         );
	f[1] = 0.5*(vol1*q1[1] + p1*dy + vol2*q2[1] + p2*dy );
	f[2] = 0.5*(vol1*q1[2] - p1*dx + vol2*q2[2] - p2*dx );
	f[3] = 0.5*(vol1*(q1[3] + p1)  + vol2*(q2[3] + p2)  );


	for (n = 0; n < 4; n++)
	{
		f[n] = f[n] + mu*(q1[n]-q2[n]);
		res1[n] = res1[n] + f[n];
		res2[n] = res2[n] - f[n];
	} 

delete[] f;

return 0;

} //end


/*********************************************************************************************************/

int flux_wall(double* x1, double* x2, double* q, double* res) 

{

//--------------------compute momentum from an individual wall face------------------------------------------
	       

int n; 
double dx, dy, ri, u,v, p;

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];

	ri = 1.0/q[0];
	u = ri*q[1];
	v = ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));
	
	res[1] = res[1] - p*dy;
	res[2] = res[2] + p*dx;
		  
        return 0;
       
} //end


/********************************************************************************************************/

int lift_wall(double* x1, double* x2, double* q, double &lift)

{


double dx, ri, u, v, p; // COMPUTE MOMENTUM FLUX FROM AN INDIVIDUAL WALL FACE

	dx = x1[0] - x2[0];
	
	ri = 1.0/q[0];
	u = ri*q[1];
	v = ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));

	lift = lift + p*dx;
	
	return 0;
       
} //end

/***************complex versions of intrinsic functions with analytic extension***************************/

double max(double arg1,double arg2)
{ 
 if(arg1 > arg2)
    return arg1 ;
 else 
    return arg2;
    
    return 0;
}

double min(double arg1, double arg2)
{
  if(arg1 <= arg2)
      return arg1;
   else
      return arg2;
   
   return 0;
}

double abs(double arg1)
{
 if(arg1 >= 0.0)
    return arg1;
  else
    return -arg1;
    
   return 0;
}
 
/*********************************************************************************************************/
