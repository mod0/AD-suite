//----------------------------------------- program airfoil------------------------------------------------
    
    
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
    
#include"const.h"
double gam, gm1, cfl, eps, mach, alpha;    
    
int main()
{
  int i;
    const int maxnode = 721801;
    const int maxcell = 720001;
    const int maxedge = 1441400;

    int **ecell = new int*[maxedge];
    for (i = 0; i < maxedge; i++)
	ecell[i] = new int[2];

    int *boun = new int[maxedge];

    int **edge = new int*[maxedge];
    for (i = 0; i < maxedge; i++)
	edge[i] =  new int[2];

    int **cell = new int*[maxcell];
    for (i = 0; i < maxcell; i++)
	cell[i] = new int[4];

    double **x = new double*[maxnode];
    for (i = 0; i < maxnode; i++)
	x[i] = new double[2];

    double **q = new double*[maxcell];
    for (i = 0; i < maxcell; i++)
	q[i] = new double[4];

    double **qold = new double*[maxcell];
    for (i = 0; i < maxcell; i++)
	qold[i] = new double[4];

    double *adt = new double[maxcell];

    double **res = new double*[maxcell];
    for (i = 0; i < maxcell; i++)
	res[i] = new double[4];

    int in1, in2, in3, in4, ic, ic1, ic2, ie, ipde, iter, niter, k, nnode, ncell, nedge;
    double lift, rms;
 
//-------------------------------read in grid and flow data--------------------------------------------
  
    input(maxnode, maxcell, maxedge, nnode, ncell, nedge, x, q, cell, edge, ecell, boun);
 
//----------------------------------main time-marching loop--------------------------------------------

    niter = 20000;

    for(iter = 0;iter < niter;++iter)
    {
 
//-----------------------------------save old flow solution--------------------------------------------
 
    for(ic=0; ic< ncell-1; ++ic)
    {
      for( ipde = 0;ipde < 4; ++ipde)
      {
        qold[ic][ipde] = q[ic][ipde];
      }
    }
 
//-------------------------------predictor/corrector update loop---------------------------------------

    for(k = 0; k < 2 ; ++k)
    { 
     for(ic = 0; ic < ncell; ++ic)
      {
        for(ipde = 0;ipde < 4; ++ipde)
         {
           res[ic][ipde] = 0.0;
         }
      }


//-----------------------------------calculate area/timstep--------------------------------------------

    for (ic = 0; ic < ncell - 1; ++ic)
    {
     in1 = cell[ic][0] -1;
     in2 = cell[ic][1] -1;
     in3 = cell[ic][2] -1;
     in4 = cell[ic][3] -1;
     time_cell(&x[in1][0], &x[in2][1], &x[in3][2],&x[in4][3],&q[ic][0],adt[ic]);
    }

    adt[ncell] = 0.0;

    //------------------------------------flux evaluation loop-------------------------------------------- 

    for (ie = 0;ie < nedge; ++ie)
    {
     in1 = edge[ie][0] -1;
     in2 = edge[ie][1] -1;
     ic1 = ecell[ie][0]-1;
     ic2 = ecell[ie][1]-1;

     if (boun[ie] == 0) 
      {
       flux_face(&x[in1][0], &x[in2][0], &q[ic1][0], &q[ic2][0],adt[ic1], adt[ic2], 
       &res[ic1][0], &res[ic2][0]);
      }
     else if (boun[ie] == 1)
      {
       flux_wall(&x[in1][0], &x[in2][0], &q[ic2][0], &res[ic2][0]);
      } 
     else if (boun[ie] == 2)
        exit(0);
    }
}

//--------------------------------------flow field update----------------------------------------------

    rms = 0.0;

    for (ic = 0;ic < ncell - 1; ++ic)
     {
      for(ipde = 0;ipde < 4; ++ipde)
       {
        q[ic][ipde] = qold[ic][ipde] - res[ic][ipde]/adt[ic];
        rms = rms + (pow((q[ic][ipde] - qold[ic][ipde]),2));
       }
     }
    }
    rms = sqrt(rms/ncell);
  

//-----------------------print iteration history, including lift calculation---------------------------

    if( (iter % 100) == 0) 
     {
      lift = 0.0;
   

    for(ie = 0;ie < nedge; ++ie)
     {
       if((boun[ie]) == 1) 
       {
        in1 = edge[ie][0] -1;
        in2 = edge[ie][1] -1;
        ic2 = ecell[ie][1]-1;
        lift_wall(&x[in1][0],&x[in2][0],&q[ic2][0], lift);
       }
     }

     printf("%d,%f,%f",iter,rms,lift);
     }


     output(ncell, q);

     for (ie = 0; ie < maxedge; ie++)
	        delete[] ecell[ie];
                delete[] ecell;

	delete[] boun;

	for (i = 0; i < maxedge; i++) 
		delete[] edge[i];
	        delete[] edge;

	for (i = 0; i < maxcell; i++)
	        delete[] cell[i];
	        delete[] cell;

	for (i = 0; i < maxnode; i++)
		delete[] x[i];
		delete[] x;

	for (i = 0; i < maxcell; i++)
		delete[] q[i];
		delete[] q;

	for (i = 0; i < maxcell; i++)
		delete[] qold[i];
		delete[] qold;

	delete[] adt;

	for (i = 0; i < maxcell; i++)
	        delete[] res[i];
		delete[] res;
return 0;
     
     
}

