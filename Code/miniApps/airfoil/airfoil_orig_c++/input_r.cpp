
//---------------------------------------INPUT.F-----------------------------------------------------------

#include <iostream>
using namespace std;

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <cstdio>

#include "const.h"

void input(int maxnode, int maxcell, int maxedge, int& nnode, int& ncell,
	   int& nedge, double **x, double** q, int** cell, int** edge, 
	   int** ecell, int* boun) 
{


        int in, ie, ic, ipde;
        double p, r, u, e, xt, yt;

//--------------------------------set universal constants-------------------------------------------------
	
        gam = 1.4; 
        gm1 = gam - 1.0;
        cfl = 0.9;
        eps = 0.05;

//--------------------------------reda in data from grid file----------------------------------------------

	ifstream f1; 
	f1.open("grid.dat");

	f1 >> nnode >> ncell >> nedge;

	ncell++;	
	if (nnode > maxnode)
		cout << "maxnode too small \n";
	if (ncell > maxcell)
		cout << "maxcell too small \n";
	if (nedge > maxedge)
		cout <<"maxedge too small \n";


	for (in = 0; in < nnode; in++)
        { 
		f1 >>  x[in][0] >> x[in][1];
	}

	for (ic = 0; ic < ncell-1; ic++) 
        {
		f1 >> cell[ic][0] >> cell[ic][1] >> cell[ic][2] >> cell[ic][3];
	}

	for (ie = 0; ie < nedge; ie++) 
        {
		f1 >>  edge[ie][0] >> edge[ie][1] >> ecell[ie][0] >> ecell[ie][1] >> boun[ie]; 
		if (boun[ie] == 2)
                {
			ecell[ie][1] = ncell;
			boun[ie] = 0;
		}	
	}

	f1.close();


//-------------------------read in data from flow file,initialising if necessary-----------------------	
	
//        ifstream f2; 
//	f2.open("flow.dat");	

//	f2 >> p >> r >> mach >> alpha;

	p = 1.0;
	r = 1.0;
        mach = 0.4;
        alpha = 3;

	alpha = alpha*atan(1.0)/45.0;

	u = sqrt(gam*p/r)*mach;
	e = p/(r*gm1) + 0.5*u*u;


	for (ic = 0; ic < ncell; ic++) 
        {
		q[ic][0] = r;
		q[ic][1] = r*u;
		q[ic][2] = 0.0;
		q[ic][3] = r*e;		
	}
/*	
        for (ic = 0; ic < ncell && !f2.eof(); ic++) 
        {
		for (ipde = 0; ipde < 4 && !f2.eof(); ipde++) 
                {
			f2 >> q[ic][ipde];	
		}
	}
*/
//	f2.close();

//-----------------------------rotate grid to specified angle of attack------------------------------------
	
        for (in = 0; in < nnode; in++) 
        { 
		xt = x[in][0];
		yt = x[in][1];
		x[in][0] = cos(alpha)*xt + sin(alpha)*yt; 
		x[in][1] = cos(alpha)*yt - sin(alpha)*xt;
	}
} 

//---------------------------------------------------------------------------------------------------------

void output(int ncell, double** q)
 {

   int ic, ipde;
   double p, r, mach, alpha;

	ifstream f3;
	f3.open("flow.dat");

	f3 >> p >> r >> mach >> alpha;

	f3.close();

	ofstream f4;
	f4.open("flowout.dat");
	f4 << " " << p << " " << r << " " << mach << " " << alpha << "\n";

	for (ic = 0; ic < ncell-1; ic++) 
        {
		for (ipde = 0; ipde < 4 ; ipde++) 
                {
		    f4 << " " << q[ic][ipde] << "\t";
		}
		f4 << "\n";
	}

	f4.close();

}
