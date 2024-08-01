/* CFD I Project 1 Code 
   Version Number: 5
   By Seyed MohammadAmin Taleghani
   Sharif University of Technology - Aerospace Engineering Department */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
using namespace std;

#define imax 101 /* Number of the Nodes in r-dir */
#define jmax 101 /* Number of the Nodes in theta-dir */
double init_guess = 0.0; /* defines the initial guess value for all T matrixes */
double err = 1e-5; /* defines the convergence critetion for iterations */

double T_jac [ imax ][jmax ], T_gauss [ imax ][jmax ], T_rsweep [ imax ][jmax ], T_thsweep [ imax ][jmax ], T_adi_rth [ imax ][jmax ];
double T_adi_thr [ imax ][jmax ], T_sor [ imax  ][jmax ], T_old [ imax ][jmax ], T_anal [ imax ][jmax ], diff [ imax ][jmax ];
double r_in = 1.0, r_out = 2.0, th_start = 0.0 , th_end = M_PI_2; /* th stands for theta */
double x [ imax ][ jmax ], y [ imax ][ jmax ], r [imax ], th [jmax ];
double dr , dth, L1, L2, L_inf, L_infty [10000], sum = 0;

int iter_jac = 0, iter_gauss = 0, iter_rsweep = 0, iter_thsweep = 0, iter_adi_rth = 0, iter_adi_thr = 0, iter_sor = 0;


/* Function Declarations */

double maxvec(double arr[], int n) 
{ 
    int i; 
      
    // Initialize maximum element 
    double max = arr[0]; 
  
    // Traverse array elements  
    // from second and compare 
    // every element with current max  
    for (i = 1; i < n; i++) 
        if (arr[i] > max) 
            max = arr[i]; 
  
    return max; 
} 


double maxarr(double arr[][jmax], int n , int m) 
{ 
    int i,j; 
      
    // Initialize maximum element 
    double max = arr[0][0]; 
  
    // Traverse array elements  
    // from second and compare 
    // every element with current max  
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
        if (arr[i][j] > max) 
			{
				max = arr[i][j];
			}
		}
	}
    return max;
}


double * thomas (double A[][imax-2], double B[])
{

	double ratio ; static double X[imax-2];

	for ( int i = 1 ; i < imax-2 ; i++  )
	{
		ratio = A[i][i-1] / A[i-1][i-1];
		A[i][i] = A[i][i] - ratio * A[i-1][i];
		B [i] = B[i] - ratio * B[i-1];
	}

	X [imax-2 -1] = B[imax-2 -1]/A[imax-2 -1][imax-2 -1];

	for ( int i = imax-2 -2 ; i >= 0 ; i--)
	{

		X [ i ] = ( B[i] - A[i][i+1]*X[i+1] )/A[i][i];

	}

	return X;
}

void initialize_all_solutions ()
{


////// Initial Guess on all nodes except boundaries////
int i, j ;

for ( i = 1; i < imax - 1; i++)
{
	for ( j = 1; j < jmax - 1; j++)
	{
			T_gauss [ i ] [ j ] = init_guess;
			T_jac   [ i ] [ j ] = init_guess;
			T_rsweep [ i ] [ j ] = init_guess;
			T_thsweep [ i ] [ j ] = init_guess;
			T_adi_rth [ i ] [ j ] = init_guess;
			T_adi_thr [ i ] [ j ] = init_guess;
			T_sor [ i ] [ j ] = init_guess;
	}
}
//

// Boundary Conditions

// @ theta = 0
for ( i = 0 ; i < imax ; i++)
{
	
	T_gauss [ i ] [ 0 ] = 0.0;
	T_jac   [ i ] [ 0 ] = 0.0;
	T_rsweep [ i ] [ 0 ] = 0.0;
	T_thsweep [ i ] [ 0 ] = 0.0;
	T_adi_rth [ i ] [ 0 ] = 0.0;
	T_adi_thr [ i ] [ 0 ] = 0.0;
	T_sor [ i ] [ 0 ] = 0.0;
}
//
// @ r = r_in
for ( j = 0 ; j < jmax ; j++)
{
	
	T_gauss [ 0 ] [ j ] = 0.0;
	T_jac   [ 0 ] [ j ] = 0.0;
	T_rsweep [ 0 ] [ j ] = 0.0;
	T_thsweep [ 0 ] [ j ] = 0.0;
	T_adi_rth [ 0 ] [ j ] = 0.0;
	T_adi_thr [ 0 ] [ j ] = 0.0;
	T_sor [ 0 ] [ j ] = 0.0;
}
//
// @ r = r_out
for ( j = 0 ; j < jmax ; j++ )
{
	
	T_gauss [ imax - 1 ] [ j ] = sin( 2 * th [ j ] );
	T_jac   [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
	T_rsweep [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
	T_thsweep [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
	T_adi_rth [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
	T_adi_thr [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
	T_sor [ imax - 1 ] [ j ] = T_gauss [ imax - 1 ] [ j ];
}
//
// @ theta = pi/2
for ( i = 0 ; i < imax ; i++ )
{	
			T_gauss [ i ] [ jmax - 1 ] = sin( M_PI*( r[i]-r_in )/( r_out-r_in ) );
			T_jac   [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
			T_rsweep [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
			T_thsweep [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
			T_adi_rth [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
			T_adi_thr [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
			T_sor [ i ] [ jmax - 1 ] = T_gauss [ i ] [ jmax - 1 ];
}
//

printf("\nAll solutions are initialized. Each node's designated temperature is: ");
cout<<init_guess<<endl<<endl;


}


void run_jacubi_method ()
{

FILE *dat_jac;
FILE *residue_jac;
FILE *norms_jac;


////// Jacobi method

//Tic///////////////////////////////////////////////////////////////////////////
std::clock_t c_start = std::clock();

int i, j;



while (1)
	{

/////// Sets T_old = T_jac//////////////////////////////////////////////////////

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_jac [ i ] [ j ];
			}
		}
///////////////////////////////////////////////////////////////////////////////

	for ( i = 1; i < imax - 1 ; i++)
	{
		for ( j = 1; j < jmax - 1 ; j++)
		{
			double term1 = (T_old[i+1][j]+T_old[i-1][j])/(dr*dr);
			double term2 = (T_old[i+1][j]-T_old[i-1][j])/(r[i]*dr*2);
			double term3 = (T_old[i][j+1]+T_old[i][j-1])/( (r[i]*dth)*(r[i]*dth) );
			double denum = 2*( 1/(dr*dr) + 1/(r[i]*dth*r[i]*dth) );

			T_jac [ i ] [ j ] =  ( term1+term2+term3 )/( denum ) ;
		}
	}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_jac [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_jac++;	

	L_infty[iter_jac] = maxarr(diff,imax,jmax);

	if (L_infty[iter_jac]<err)
		{

			break;
		}
	
	}

//Toc //////////////////////////////////////////////////////////////////////
 std::clock_t c_end = std::clock();
 long double cpu_time_jac = 1000 * (c_end-c_start) / CLOCKS_PER_SEC;


//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_jac [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_jac [ i ] [ j ] - T_anal [ i ] [ j ]) );
				L1 = sum / (imax*jmax);
			}
		}
	
sum = 0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_jac [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
				L2 = pow ( sum / (imax*jmax) , 0.5);
			}
		}


/* Tecplot Temperatures Outputs for Jacubi*/
/* Header */
dat_jac=fopen("dat_jac.dat","w");
fprintf ( dat_jac,"Variables = " );
fprintf ( dat_jac, "x");
fprintf ( dat_jac, "," );
fprintf ( dat_jac, "y");
fprintf ( dat_jac, "," );
fprintf ( dat_jac, "T");
fprintf ( dat_jac, "\n");
fprintf ( dat_jac, "Zone I = " );
fprintf ( dat_jac, "%d", imax );
fprintf ( dat_jac, ", J = " );
fprintf ( dat_jac, "%d", jmax );
fprintf ( dat_jac, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_jac, "%e", x [i][j]);
		fprintf (dat_jac, " " );
		fprintf (dat_jac, "%e", y [i][j]);
		fprintf (dat_jac, " " );
		fprintf (dat_jac, "%e", T_jac [i][j]);
		fprintf (dat_jac, "\n");
	}
}
fclose(dat_jac);

/* Tecplot Residue Outputs for Jacubi*/
/* Header */
residue_jac=fopen("residue_jac.dat","w");
fprintf ( residue_jac,"Variables = " );
fprintf ( residue_jac, "i");
fprintf ( residue_jac, "," );
fprintf ( residue_jac, "e");
fprintf ( residue_jac, "\n");
fprintf ( residue_jac, "Zone I = " );
fprintf ( residue_jac, "%d", iter_jac );
fprintf ( residue_jac, "\n" );
/* Data */
for ( i = 1; i <= iter_jac; i++ )
	{
		fprintf (residue_jac, "%d", i);
		fprintf (residue_jac, " " );
		fprintf (residue_jac, "%e", log10(L_infty[i]) );
		fprintf (residue_jac, "\n");
	}
fclose(residue_jac);

/* Error Norms Outputs */
norms_jac=fopen("norms_jac.txt","w");
fprintf ( norms_jac, "Number of Iteration for the Jacubi Method = " );
fprintf ( norms_jac,"%d", iter_jac );
fprintf ( norms_jac, "\n" );
fprintf ( norms_jac, "CPU Time (ms) = " );
fprintf ( norms_jac, "%Lf", cpu_time_jac );
fprintf ( norms_jac, "\n" );
fprintf ( norms_jac, "L_infty = " );
fprintf ( norms_jac, "%e", L_inf );
fprintf ( norms_jac, "\n");
fprintf ( norms_jac, "L_1 = " );
fprintf ( norms_jac, "%e", L1 );
fprintf ( norms_jac, "\n");
fprintf ( norms_jac, "L_2 = " );
fprintf ( norms_jac, "%e", L2 );
fclose(norms_jac);


printf("Numerical solution using Jacubi iteration has finished. Number of iterations for this method: %d\n\n",iter_jac);

}




void run_gauss_method ()
{

FILE *dat_gauss;
FILE *residue_gauss;
FILE *norms_gauss;


////// Gauss-Seidel method
// Tic
std::clock_t c_start = std::clock();

int i,j;

while (1)

	{

/////// Sets T_old = T_gauss

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_gauss [ i ] [ j ];
			}
		}
	
///////  main loop starts here
	for ( i = 1; i < imax - 1 ; i++)
	{
		for ( j = 1; j < jmax - 1 ; j++)
		{
			double term1 = (T_old[i+1][j]+T_gauss[i-1][j])/(dr*dr);
			double term2 = (T_old[i+1][j]-T_gauss[i-1][j])/(r[i]*dr*2);
			double term3 = (T_old[i][j+1]+T_gauss[i][j-1])/( (r[i]*dth)*(r[i]*dth) );
			double denum = 2.0*( 1/(dr*dr) + 1/(r[i]*dth*r[i]*dth) );

			T_gauss [ i ] [ j ] =  ( term1+term2+term3 )/( denum );
		}
	}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_gauss [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_gauss++;

	L_infty[iter_gauss] = maxarr(diff,imax,jmax);

	if (L_infty[iter_gauss]<err)
		{

			break;
		}
	
	}


//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_gauss [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_gauss [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);

	
sum = 0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_gauss [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


// Toc
std::clock_t c_end = std::clock();
long double cpu_time_gauss = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

/* Tecplot Temperatures Outputs for Gauss-Seidle */
/* Header */
dat_gauss=fopen("dat_gauss.dat","w");
fprintf ( dat_gauss,"Variables = " );
fprintf ( dat_gauss, "x");
fprintf ( dat_gauss, "," );
fprintf ( dat_gauss, "y");
fprintf ( dat_gauss, "," );
fprintf ( dat_gauss, "T");
fprintf ( dat_gauss, "\n");
fprintf ( dat_gauss, "Zone I = " );
fprintf ( dat_gauss, "%d", imax );
fprintf ( dat_gauss, ", J = " );
fprintf ( dat_gauss, "%d", jmax );
fprintf ( dat_gauss, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_gauss, "%e", x [i][j]);
		fprintf (dat_gauss, " " );
		fprintf (dat_gauss, "%e", y [i][j]);
		fprintf (dat_gauss, " " );
		fprintf (dat_gauss, "%e", T_gauss [i][j]);
		fprintf (dat_gauss, "\n");
	}
}
fclose(dat_gauss);

/* Tecplot Residue Outputs for Gauss-Seidle*/
/* Header */
residue_gauss=fopen("residue_gauss.dat","w");
fprintf ( residue_gauss,"Variables = " );
fprintf ( residue_gauss, "i");
fprintf ( residue_gauss, "," );
fprintf ( residue_gauss, "e");
fprintf ( residue_gauss, "\n");
fprintf ( residue_gauss, "Zone I = " );
fprintf ( residue_gauss, "%d", iter_gauss );
fprintf ( residue_gauss, "\n" );
/* Data */
for ( i = 1; i <= iter_gauss; i++ )
	{
		fprintf (residue_gauss, "%d", i);
		fprintf (residue_gauss, " " );
		fprintf (residue_gauss, "%e", log10(L_infty[i]) );
		fprintf (residue_gauss, "\n");
	}
fclose(residue_gauss);

/* Error Norms Outputs */
norms_gauss=fopen("norms_gauss.txt","w");
fprintf ( norms_gauss, "Number of Iteration for the Gauss-Seidle Method = " );
fprintf ( norms_gauss,"%d", iter_gauss );
fprintf ( norms_gauss, "\n" );
fprintf ( norms_gauss, "CPU Time (ms) = " );
fprintf ( norms_gauss, "%Lf", cpu_time_gauss );
fprintf ( norms_gauss, "\n" );
fprintf ( norms_gauss, "L_infty = " );
fprintf ( norms_gauss, "%e", L_inf );
fprintf ( norms_gauss, "\n");
fprintf ( norms_gauss, "L_1 = " );
fprintf ( norms_gauss, "%e", L1 );
fprintf ( norms_gauss, "\n");
fprintf ( norms_gauss, "L_2 = " );
fprintf ( norms_gauss, "%e", L2 );
fclose(norms_gauss);


printf("Numerical solution using Gauss-Seidle iteration has finished. Number of iterations for this method: %d\n\n",iter_gauss);

}


void run_rsweep_method ()
{
	FILE *dat_rsweep;
	FILE *residue_rsweep;
	FILE *norms_rsweep;

////////// r-sweep Method
////////// Tic
	std::clock_t c_start = std::clock();

	int i,j;

//////// Defines A, B and X matrices for A*X=B;
	double A [imax-2][imax-2] , B [imax-2];
	double	*X;

	for ( i = 0; i < imax - 2; i++)
	{
	for ( j = 0; j < jmax - 2; j++)
		{
			A [ i ] [ j ] = 0.0;
		}
	}


///////  main loop starts here

	while (1)
	{

////////////// Sets T_old = T_rsweep

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_rsweep [ i ] [ j ];
			}
		}

for ( j = 1 ; j < jmax-1 ; j++ )
		{

//////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2.0 * ( 1/(dr*dr) + 1/(pow(r[1]*dth,2)) );
			A [0][1] = 1.0/(dr*dr) + 1/(2.0*r[1]*dr) ;
			B [0]    = - ( T_rsweep [1][j-1] + T_rsweep [1][j+1] )/(pow(r[1]*dth,2)) - ( 1.0/(dr*dr) - 1.0/(2.0*r[1]*dr) )*T_rsweep[0][j];

			for ( i = 1 ; i < imax-2 -1 ; i++)
			{

				A [i][i-1] = 1.0/(dr*dr) - 1.0/(2*r[i+1]*dr);
				A [i][i]   = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[i+1]*dth,2)) );
				A [i][i+1] = 1.0/(dr*dr) + 1/(2*r[i+1]*dr);
				B [i]      = - ( T_rsweep [i+1][j-1] + T_rsweep [i+1][j+1] )/(pow(r[i+1]*dth,2));

			}

			A [imax-2 -1][imax-2 -2] = 1.0/(dr*dr) - 1.0/(2*r[imax-2]*dr);
			A [imax-2 -1][imax-2 -1] = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[imax-2]*dth,2)) );
			B [imax-2 -1]			 = - ( T_rsweep [imax-2][j-1] + T_rsweep [imax-2][j+1] )/(pow(r[imax-2]*dth,2)) - (1.0/(dr*dr) + 1.0/(2.0*r[imax-2]*dr))*T_rsweep[imax-1][j];

			X = thomas (A,B);


			for ( i = 0 ; i < imax-2 ; i++)
			{
				T_rsweep [i+1][j] = *(X+i);
			}

		}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_rsweep [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_rsweep++;	

	L_infty[iter_rsweep] = maxarr(diff,imax,jmax);

	if (L_infty[iter_rsweep]<err)
		{

			break;
		}

	}

	// Toc
	std::clock_t c_end = std::clock();
	long double cpu_time_rsweep = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_rsweep [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_rsweep [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);
	

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_rsweep [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


/* Tecplot Temperatures Outputs for r-sweep method */
/* Header */
dat_rsweep=fopen("dat_rsweep.dat","w");
fprintf ( dat_rsweep,"Variables = " );
fprintf ( dat_rsweep, "x");
fprintf ( dat_rsweep, "," );
fprintf ( dat_rsweep, "y");
fprintf ( dat_rsweep, "," );
fprintf ( dat_rsweep, "T");
fprintf ( dat_rsweep, "\n");
fprintf ( dat_rsweep, "Zone I = " );
fprintf ( dat_rsweep, "%d", imax );
fprintf ( dat_rsweep, ", J = " );
fprintf ( dat_rsweep, "%d", jmax );
fprintf ( dat_rsweep, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_rsweep, "%e", x [i][j]);
		fprintf (dat_rsweep, " " );
		fprintf (dat_rsweep, "%e", y [i][j]);
		fprintf (dat_rsweep, " " );
		fprintf (dat_rsweep, "%e", T_rsweep [i][j]);
		fprintf (dat_rsweep, "\n");
	}
}
fclose(dat_rsweep);

/* Tecplot Residue Outputs for r-sweep Method */
/* Header */
residue_rsweep=fopen("residue_rsweep.dat","w");
fprintf ( residue_rsweep,"Variables = " );
fprintf ( residue_rsweep, "i");
fprintf ( residue_rsweep, "," );
fprintf ( residue_rsweep, "e");
fprintf ( residue_rsweep, "\n");
fprintf ( residue_rsweep, "Zone I = " );
fprintf ( residue_rsweep, "%d", iter_rsweep );
fprintf ( residue_rsweep, "\n" );
/* Data */
for ( i = 1; i <= iter_rsweep; i++ )
	{
		fprintf (residue_rsweep, "%d", i);
		fprintf (residue_rsweep, " " );
		fprintf (residue_rsweep, "%e", log10(L_infty[i]) );
		fprintf (residue_rsweep, "\n");
	}
fclose(residue_rsweep);

/* Error Norms Outputs */
norms_rsweep=fopen("norms_rsweep.txt","w");
fprintf ( norms_rsweep, "Number of Iteration for the r-sweep Method = " );
fprintf ( norms_rsweep,"%d", iter_rsweep );
fprintf ( norms_rsweep, "\n" );
fprintf ( norms_rsweep, "CPU Time (ms) = " );
fprintf ( norms_rsweep, "%Lf", cpu_time_rsweep );
fprintf ( norms_rsweep, "\n" );
fprintf ( norms_rsweep, "L_infty = " );
fprintf ( norms_rsweep, "%e", L_inf );
fprintf ( norms_rsweep, "\n");
fprintf ( norms_rsweep, "L_1 = " );
fprintf ( norms_rsweep, "%e", L1 );
fprintf ( norms_rsweep, "\n");
fprintf ( norms_rsweep, "L_2 = " );
fprintf ( norms_rsweep, "%e", L2 );
fclose(norms_rsweep);


printf("Numerical solution using r-sweep block iteration has finished. Number of iterations for this method: %d\n\n",iter_rsweep);


}


void run_thsweep_method ()
{
	FILE *dat_thsweep;
	FILE *residue_thsweep;
	FILE *norms_thsweep;


////////// theta-sweep Method
////////// Tic
	std::clock_t c_start = std::clock();

	int i,j;

//////// Defines A, B and X matrices for A*X=B;
	double A [jmax-2][jmax-2] , B [jmax-2];
	double	*X;

	for ( i = 0; i < jmax - 2; i++)
	{
	for ( j = 0; j < jmax - 2; j++)
		{
			A [ i ] [ j ] = 0.0;
		}
	}


///////  main loop starts here

	while (1)
	{

////////////// Sets T_old = T_thsweep

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_thsweep [ i ] [ j ];
			}
	}

///////////////////////////////////////////

		for ( i = 1 ; i < imax -1 ; i++ )
		{

//////////////// Fills A, B and X matrices (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			A [0][1] = 1/(pow(r[i]*dth,2)) ;
			B [0]    = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_thsweep [i-1][1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_thsweep [i+1][1]  - (1/(pow(r[i]*dth,2)))*T_thsweep[i][0];

			for ( j = 1 ; j < jmax -2 -1; j++)
			{

				A [j][j-1] = 1/(pow(r[i]*dth,2));
				A [j][j]   = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
				A [j][j+1] = 1/(pow(r[i]*dth,2));
				B [j]      = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_thsweep [i-1][j+1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_thsweep [i+1][j+1];

			}

			A [imax-2 -1][imax-2 -2] = 1/(pow(r[i]*dth,2));
			A [imax-2 -1][imax-2 -1] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			B [imax-2 -1]			 = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_thsweep [i-1][jmax-2] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_thsweep [i+1][jmax-2]  - (1/(pow(r[i]*dth,2)))*T_thsweep[i][jmax-1];

			X = thomas (A,B);


			for ( j = 0 ; j < jmax-2 ; j++)
			{

				T_thsweep [i][j+1] = *(X + j);
			}

		}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_thsweep [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_thsweep++;	

	L_infty[iter_thsweep] = maxarr(diff,imax,jmax);

	if (L_infty[iter_thsweep]<err)
		{

			break;
		}

	}

	// Toc
	std::clock_t c_end = std::clock();
	long double cpu_time_thsweep = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_thsweep [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_thsweep [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);
	
sum = 0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_thsweep [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


/* Tecplot Temperatures Outputs for theta-sweep */
/* Header */
dat_thsweep=fopen("dat_thsweep.dat","w");
fprintf ( dat_thsweep,"Variables = " );
fprintf ( dat_thsweep, "x");
fprintf ( dat_thsweep, "," );
fprintf ( dat_thsweep, "y");
fprintf ( dat_thsweep, "," );
fprintf ( dat_thsweep, "T");
fprintf ( dat_thsweep, "\n");
fprintf ( dat_thsweep, "Zone I = " );
fprintf ( dat_thsweep, "%d", imax );
fprintf ( dat_thsweep, ", J = " );
fprintf ( dat_thsweep, "%d", jmax );
fprintf ( dat_thsweep, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_thsweep, "%e", x [i][j]);
		fprintf (dat_thsweep, " " );
		fprintf (dat_thsweep, "%e", y [i][j]);
		fprintf (dat_thsweep, " " );
		fprintf (dat_thsweep, "%e", T_thsweep [i][j]);
		fprintf (dat_thsweep, "\n");
	}
}
fclose(dat_thsweep);

/* Tecplot Residue Outputs for theta-sweep*/
/* Header */
residue_thsweep=fopen("residue_thsweep.dat","w");
fprintf ( residue_thsweep,"Variables = " );
fprintf ( residue_thsweep, "i");
fprintf ( residue_thsweep, "," );
fprintf ( residue_thsweep, "e");
fprintf ( residue_thsweep, "\n");
fprintf ( residue_thsweep, "Zone I = " );
fprintf ( residue_thsweep, "%d", iter_thsweep );
fprintf ( residue_thsweep, "\n" );
/* Data */
for ( i = 1; i <= iter_thsweep; i++ )
	{
		fprintf (residue_thsweep, "%d", i);
		fprintf (residue_thsweep, " " );
		fprintf (residue_thsweep, "%e", log10(L_infty[i]) );
		fprintf (residue_thsweep, "\n");
	}
fclose(residue_thsweep);

/* Error Norms Outputs */
norms_thsweep=fopen("norms_thsweep.txt","w");
fprintf ( norms_thsweep, "Number of Iteration for the theta-sweep Method = " );
fprintf ( norms_thsweep,"%d", iter_thsweep );
fprintf ( norms_thsweep, "\n" );
fprintf ( norms_thsweep, "CPU Time (ms) = " );
fprintf ( norms_thsweep, "%Lf", cpu_time_thsweep );
fprintf ( norms_thsweep, "\n" );
fprintf ( norms_thsweep, "L_infty = " );
fprintf ( norms_thsweep, "%e", L_inf );
fprintf ( norms_thsweep, "\n");
fprintf ( norms_thsweep, "L_1 = " );
fprintf ( norms_thsweep, "%e", L1 );
fprintf ( norms_thsweep, "\n");
fprintf ( norms_thsweep, "L_2 = " );
fprintf ( norms_thsweep, "%e", L2 );
fclose(norms_thsweep);


printf("Numerical solution using theta-sweep block iteration has finished. Number of iterations for this method: %d\n\n",iter_thsweep);

}


void run_adi_rth_method ()

{
	FILE *daT_adi_rth;
	FILE *residue_adi_rth;
	FILE *norms_adi_rth;


/////////// ADI Method
////////// Tic
	std::clock_t c_start = std::clock();

	int i,j;

//////// Defines A, B and X matrices for A*X=B;
	double A [imax-2][imax-2] , B [imax-2];
	double	*X;

	for ( i = 0; i < imax - 2; i++)
	{
	for ( j = 0; j < jmax - 2; j++)
		{
			A [ i ] [ j ] = 0.0;
		}
	}


///////  main loop starts here

	while (1)
	{

////////////// Sets T_old = T_adi_rth

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_adi_rth [ i ] [ j ];
			}
		}

/////////////////  r-sweep loop   ////////////////////////////

		for ( j = 1 ; j < jmax-1 ; j++ )
		{

//////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2.0 * ( 1/(dr*dr) + 1/(pow(r[1]*dth,2)) );
			A [0][1] = 1.0/(dr*dr) + 1/(2.0*r[1]*dr) ;
			B [0]    = - ( T_adi_rth [1][j-1] + T_adi_rth [1][j+1] )/(pow(r[1]*dth,2)) - ( 1.0/(dr*dr) - 1.0/(2.0*r[1]*dr) )*T_adi_rth[0][j];

			for ( i = 1 ; i < imax-2 -1 ; i++)
			{

				A [i][i-1] = 1.0/(dr*dr) - 1.0/(2*r[i+1]*dr);
				A [i][i]   = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[i+1]*dth,2)) );
				A [i][i+1] = 1.0/(dr*dr) + 1/(2*r[i+1]*dr);
				B [i]      = - ( T_adi_rth [i+1][j-1] + T_adi_rth [i+1][j+1] )/(pow(r[i+1]*dth,2));

			}

			A [imax-2 -1][imax-2 -2] = 1.0/(dr*dr) - 1.0/(2*r[imax-2]*dr);
			A [imax-2 -1][imax-2 -1] = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[imax-2]*dth,2)) );
			B [imax-2 -1]			 = - ( T_adi_rth [imax-2][j-1] + T_adi_rth [imax-2][j+1] )/(pow(r[imax-2]*dth,2)) - (1.0/(dr*dr) + 1.0/(2.0*r[imax-2]*dr))*T_adi_rth[imax-1][j];

			X = thomas (A,B);


			for ( i = 0 ; i < imax-2 ; i++)
			{
				T_adi_rth [i+1][j] = *(X+i);
			}

		}
/////////////////  r-sweep loop ends   ////////////////////////////



/////////////////  theta-sweep loop    ////////////////////////////
		for ( i = 1 ; i < imax -1 ; i++ )
		{

//////////////// Fills A, B and X matrices for theta-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			A [0][1] = 1/(pow(r[i]*dth,2)) ;
			B [0]    = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_rth [i-1][1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_rth [i+1][1]  - (1/(pow(r[i]*dth,2)))*T_adi_rth[i][0];

			for ( j = 1 ; j < jmax -2 -1; j++)
			{

				A [j][j-1] = 1/(pow(r[i]*dth,2));
				A [j][j]   = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
				A [j][j+1] = 1/(pow(r[i]*dth,2));
				B [j]      = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_rth [i-1][j+1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_rth [i+1][j+1];

			}

			A [imax-2 -1][imax-2 -2] = 1/(pow(r[i]*dth,2));
			A [imax-2 -1][imax-2 -1] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			B [imax-2 -1]			 = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_rth [i-1][jmax-2] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_rth [i+1][jmax-2]  - (1/(pow(r[i]*dth,2)))*T_adi_rth[i][jmax-1];

			X = thomas (A,B);


			for ( j = 0 ; j < jmax-2 ; j++)
			{

				T_adi_rth [i][j+1] = *(X + j);
			}

/////////////////  theta-sweep loop ends   ////////////////////////////

		}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_adi_rth [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_adi_rth++;	

	L_infty[iter_adi_rth] = maxarr(diff,imax,jmax);

	if (L_infty[iter_adi_rth]<err)
		{

			break;
		}

}

// Toc
	std::clock_t c_end = std::clock();
	long double cpu_time_adi = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;


//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_adi_rth [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_adi_rth [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);
	

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_adi_rth [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


/* Tecplot Temperatures Outputs for r-sweep method */
/* Header */
daT_adi_rth=fopen("daT_adi_rth.dat","w");
fprintf ( daT_adi_rth,"Variables = " );
fprintf ( daT_adi_rth, "x");
fprintf ( daT_adi_rth, "," );
fprintf ( daT_adi_rth, "y");
fprintf ( daT_adi_rth, "," );
fprintf ( daT_adi_rth, "T");
fprintf ( daT_adi_rth, "\n");
fprintf ( daT_adi_rth, "Zone I = " );
fprintf ( daT_adi_rth, "%d", imax );
fprintf ( daT_adi_rth, ", J = " );
fprintf ( daT_adi_rth, "%d", jmax );
fprintf ( daT_adi_rth, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (daT_adi_rth, "%e", x [i][j]);
		fprintf (daT_adi_rth, " " );
		fprintf (daT_adi_rth, "%e", y [i][j]);
		fprintf (daT_adi_rth, " " );
		fprintf (daT_adi_rth, "%e", T_adi_rth [i][j]);
		fprintf (daT_adi_rth, "\n");
	}
}
fclose(daT_adi_rth);

/* Tecplot Residue Outputs for ADI method*/
/* Header */
residue_adi_rth=fopen("residue_adi_rth.dat","w");
fprintf ( residue_adi_rth,"Variables = " );
fprintf ( residue_adi_rth, "i");
fprintf ( residue_adi_rth, "," );
fprintf ( residue_adi_rth, "e");
fprintf ( residue_adi_rth, "\n");
fprintf ( residue_adi_rth, "Zone I = " );
fprintf ( residue_adi_rth, "%d", iter_adi_rth );
fprintf ( residue_adi_rth, "\n" );
/* Data */
for ( i = 1; i <= iter_adi_rth; i++ )
	{
		fprintf (residue_adi_rth, "%d", i);
		fprintf (residue_adi_rth, " " );
		fprintf (residue_adi_rth, "%e", log10(L_infty[i]) );
		fprintf (residue_adi_rth, "\n");
	}
fclose(residue_adi_rth);

/* Error Norms Outputs */
norms_adi_rth=fopen("norms_adi_rth.txt","w");
fprintf ( norms_adi_rth, "Number of Iteration for the ADI (r-sweep then theta sweep) Method = " );
fprintf ( norms_adi_rth,"%d", iter_adi_rth );
fprintf ( norms_adi_rth, "\n" );
fprintf ( norms_adi_rth, "CPU Time (ms) = " );
fprintf ( norms_adi_rth, "%Lf", cpu_time_adi );
fprintf ( norms_adi_rth, "\n" );
fprintf ( norms_adi_rth, "L_infty = " );
fprintf ( norms_adi_rth, "%e", L_inf );
fprintf ( norms_adi_rth, "\n");
fprintf ( norms_adi_rth, "L_1 = " );
fprintf ( norms_adi_rth, "%e", L1 );
fprintf ( norms_adi_rth, "\n");
fprintf ( norms_adi_rth, "L_2 = " );
fprintf ( norms_adi_rth, "%e", L2 );
fclose(norms_adi_rth);


printf("Numerical solution using ADI block (r-sweep then theta sweep) iteration has finished. Number of iterations for this method: %d\n\n",iter_adi_rth);


}

void run_adi_thr_method ()

{
	FILE *dat_adi_thr;
	FILE *residue_adi_thr;
	FILE *norms_adi_thr;

/////////// ADI Method
////////// Tic
	std::clock_t c_start = std::clock();

	int i,j;

//////// Defines A, B and X matrices for A*X=B;
	double A [imax-2][imax-2] , B [imax-2];
	double	*X;

	for ( i = 0; i < imax - 2; i++)
	{
	for ( j = 0; j < jmax - 2; j++)
		{
			A [ i ] [ j ] = 0.0;
		}
	}


///////  main loop starts here

	while (1)
	{

////////////// Sets T_old = T_adi_rth

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_adi_thr [ i ] [ j ];
			}
		}


/////////////////  theta-sweep loop    ////////////////////////////
		for ( i = 1 ; i < imax -1 ; i++ )
		{

//////////////// Fills A, B and X matrices for theta-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			A [0][1] = 1/(pow(r[i]*dth,2)) ;
			B [0]    = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_thr [i-1][1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_thr [i+1][1]  - (1/(pow(r[i]*dth,2)))*T_adi_thr[i][0];

			for ( j = 1 ; j < jmax -2 -1; j++)
			{

				A [j][j-1] = 1/(pow(r[i]*dth,2));
				A [j][j]   = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
				A [j][j+1] = 1/(pow(r[i]*dth,2));
				B [j]      = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_thr [i-1][j+1] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_thr [i+1][j+1];

			}

			A [imax-2 -1][imax-2 -2] = 1/(pow(r[i]*dth,2));
			A [imax-2 -1][imax-2 -1] = -2 * ( 1/(dr*dr) + 1/(pow(r[i]*dth,2)) );
			B [imax-2 -1]			 = - ( 1/(dr*dr) - 1/(2*r[i]*dr) )*T_adi_thr [i-1][jmax-2] - ( 1/(dr*dr) + 1/(2*r[i]*dr) )*T_adi_thr [i+1][jmax-2]  - (1/(pow(r[i]*dth,2)))*T_adi_thr[i][jmax-1];

			X = thomas (A,B);


			for ( j = 0 ; j < jmax-2 ; j++)
			{

				T_adi_thr [i][j+1] = *(X + j);
			}

		}
/////////////////  theta-sweep loop ends   ////////////////////////////


/////////////////////   r-sweep loop   ////////////////////////////////

		for ( j = 1 ; j < jmax-1 ; j++ )
		{

//////////// Fills A, B and X matrices for r-sweep (the 1st and last row formulas are different and are filled outside the loop);

			A [0][0] = -2.0 * ( 1/(dr*dr) + 1/(pow(r[1]*dth,2)) );
			A [0][1] = 1.0/(dr*dr) + 1/(2.0*r[1]*dr) ;
			B [0]    = - ( T_adi_thr [1][j-1] + T_adi_thr [1][j+1] )/(pow(r[1]*dth,2)) - ( 1.0/(dr*dr) - 1.0/(2.0*r[1]*dr) )*T_adi_thr[0][j];

			for ( i = 1 ; i < imax-2 -1 ; i++)
			{

				A [i][i-1] = 1.0/(dr*dr) - 1.0/(2*r[i+1]*dr);
				A [i][i]   = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[i+1]*dth,2)) );
				A [i][i+1] = 1.0/(dr*dr) + 1/(2*r[i+1]*dr);
				B [i]      = - ( T_adi_thr [i+1][j-1] + T_adi_thr [i+1][j+1] )/(pow(r[i+1]*dth,2));

			}

			A [imax-2 -1][imax-2 -2] = 1.0/(dr*dr) - 1.0/(2*r[imax-2]*dr);
			A [imax-2 -1][imax-2 -1] = -2.0 * ( 1.0/(dr*dr) + 1.0/(pow(r[imax-2]*dth,2)) );
			B [imax-2 -1]			 = - ( T_adi_thr [imax-2][j-1] + T_adi_thr [imax-2][j+1] )/(pow(r[imax-2]*dth,2)) - (1.0/(dr*dr) + 1.0/(2.0*r[imax-2]*dr))*T_adi_thr[imax-1][j];

			X = thomas (A,B);


			for ( i = 0 ; i < imax-2 ; i++)
			{
				T_adi_thr [i+1][j] = *(X+i);
			}

		}

/////////////////  r-sweep loop ends   ////////////////////////////


	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_adi_thr [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_adi_thr++;	

	L_infty[iter_adi_thr] = maxarr(diff,imax,jmax);

	if (L_infty[iter_adi_thr]<err)
		{

			break;
		}

}

// Toc
	std::clock_t c_end = std::clock();
	long double cpu_time_adi = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

//Calculates L1 and L2 and L_infty from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_adi_thr [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_adi_thr [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);
	

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_adi_thr [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


/* Tecplot Temperatures Outputs for r-sweep method */
/* Header */
dat_adi_thr=fopen("dat_adi_rth.dat","w");
fprintf ( dat_adi_thr,"Variables = " );
fprintf ( dat_adi_thr, "x");
fprintf ( dat_adi_thr, "," );
fprintf ( dat_adi_thr, "y");
fprintf ( dat_adi_thr, "," );
fprintf ( dat_adi_thr, "T");
fprintf ( dat_adi_thr, "\n");
fprintf ( dat_adi_thr, "Zone I = " );
fprintf ( dat_adi_thr, "%d", imax );
fprintf ( dat_adi_thr, ", J = " );
fprintf ( dat_adi_thr, "%d", jmax );
fprintf ( dat_adi_thr, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_adi_thr, "%e", x [i][j]);
		fprintf (dat_adi_thr, " " );
		fprintf (dat_adi_thr, "%e", y [i][j]);
		fprintf (dat_adi_thr, " " );
		fprintf (dat_adi_thr, "%e", T_adi_thr [i][j]);
		fprintf (dat_adi_thr, "\n");
	}
}
fclose(dat_adi_thr);

/* Tecplot Residue Outputs for ADI method*/
/* Header */
residue_adi_thr=fopen("residue_adi_thr.dat","w");
fprintf ( residue_adi_thr,"Variables = " );
fprintf ( residue_adi_thr, "i");
fprintf ( residue_adi_thr, "," );
fprintf ( residue_adi_thr, "e");
fprintf ( residue_adi_thr, "\n");
fprintf ( residue_adi_thr, "Zone I = " );
fprintf ( residue_adi_thr, "%d", iter_adi_thr );
fprintf ( residue_adi_thr, "\n" );
/* Data */
for ( i = 1; i <= iter_adi_rth; i++ )
	{
		fprintf (residue_adi_thr, "%d", i);
		fprintf (residue_adi_thr, " " );
		fprintf (residue_adi_thr, "%e", log10(L_infty[i]) );
		fprintf (residue_adi_thr, "\n");
	}
fclose(residue_adi_thr);

/* Error Norms Outputs */
norms_adi_thr=fopen("norms_adi_thr.txt","w");
fprintf ( norms_adi_thr, "Number of Iteration for the ADI (theta sweep then r-sweep) Method = " );
fprintf ( norms_adi_thr,"%d", iter_adi_thr );
fprintf ( norms_adi_thr, "\n" );
fprintf ( norms_adi_thr, "CPU Time (ms) = " );
fprintf ( norms_adi_thr, "%Lf", cpu_time_adi );
fprintf ( norms_adi_thr, "\n" );
fprintf ( norms_adi_thr, "L_infty = " );
fprintf ( norms_adi_thr, "%e", L_inf );
fprintf ( norms_adi_thr, "\n");
fprintf ( norms_adi_thr, "L_1 = " );
fprintf ( norms_adi_thr, "%e", L1 );
fprintf ( norms_adi_thr, "\n");
fprintf ( norms_adi_thr, "L_2 = " );
fprintf ( norms_adi_thr, "%e", L2 );
fclose(norms_adi_thr);


printf("Numerical solution using ADI block (theta sweep then r-sweep) iteration has finished. Number of iterations for this method: %d\n\n",iter_adi_thr);


}

void run_sor_method ()
{
	FILE *dat_sor;
	FILE *residue_sor;
	FILE *norms_sor;


	double omega = 2*(1-(M_PI/imax));

////////// Point Gauss-Seidle SOR method
// Tic
std::clock_t c_start = std::clock();

int i,j;

while (1)

	{

/////// Sets T_old = T_sor

	for ( i = 0; i < imax ; i++)
	{	
		for ( j = 0; j < jmax ; j++)
			{
				T_old [ i ] [ j ] = T_sor [ i ] [ j ];
			}
		}
	
///////  main loop starts here
	for ( i = 1; i < imax - 1 ; i++)
	{
		for ( j = 1; j < jmax - 1 ; j++)
		{

			double term1 = (T_old[i+1][j]+T_sor[i-1][j])/(dr*dr);
			double term2 = (T_old[i+1][j]-T_sor[i-1][j])/(r[i]*dr*2);
			double term3 = (T_old[i][j+1]+T_sor[i][j-1])/( (r[i]*dth)*(r[i]*dth) );
			double coeff = 1/( (2.0*( 1/(dr*dr) + 1/(r[i]*dth*r[i]*dth))) );

			T_sor [ i ] [ j ] =  (1.0-omega)*T_old[i][j] + omega*coeff*( term1+term2+term3 );
		}
	}

	for ( i = 0; i < imax  ; i++)
	{
		for ( j = 0; j < jmax  ; j++)
		{
			diff   [ i ] [ j ] =  abs ( (T_sor [ i ] [ j ] - T_old [ i ] [ j ]) );
		}
	
	}

	iter_sor++;

	L_infty[iter_sor] = maxarr(diff,imax,jmax);

	if (L_infty[iter_sor]<err)
		{

			break;
		}
	
	}

//Calculates L1 and L2 and L_inf from the last iteration //////////////


for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				diff   [ i ] [ j ] =  abs ( (T_sor [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}

L_inf = maxarr(diff,imax,jmax);

sum = 0.0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + abs ( (T_sor [ i ] [ j ] - T_anal [ i ] [ j ]) );
			}
		}
L1 = sum / (imax*jmax);

	
sum = 0;
for ( i = 0 ; i<imax ; i++)
		{
	for ( j = 0 ; j<jmax; j++)
			{
				sum = sum + pow ( (T_sor [ i ] [ j ] - T_anal [ i ] [ j ]),2 );
			}
		}
L2 = pow ( sum / (imax*jmax) , 0.5);


// Toc
std::clock_t c_end = std::clock();
long double cpu_time_sor = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

/* Tecplot Temperatures Outputs for Point Gauss-Seidle SOR */
/* Header */
dat_sor=fopen("dat_sor.dat","w");
fprintf ( dat_sor,"Variables = " );
fprintf ( dat_sor, "x");
fprintf ( dat_sor, "," );
fprintf ( dat_sor, "y");
fprintf ( dat_sor, "," );
fprintf ( dat_sor, "T");
fprintf ( dat_sor, "\n");
fprintf ( dat_sor, "Zone I = " );
fprintf ( dat_sor, "%d", imax );
fprintf ( dat_sor, ", J = " );
fprintf ( dat_sor, "%d", jmax );
fprintf ( dat_sor, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_sor, "%e", x [i][j]);
		fprintf (dat_sor, " " );
		fprintf (dat_sor, "%e", y [i][j]);
		fprintf (dat_sor, " " );
		fprintf (dat_sor, "%e", T_sor [i][j]);
		fprintf (dat_sor, "\n");
	}
}
fclose(dat_sor);

/* Tecplot Residue Outputs for Point Gauss-Seidle SOR*/
/* Header */
residue_sor=fopen("residue_sor.dat","w");
fprintf ( residue_sor,"Variables = " );
fprintf ( residue_sor, "i");
fprintf ( residue_sor, "," );
fprintf ( residue_sor, "e");
fprintf ( residue_sor, "\n");
fprintf ( residue_sor, "Zone I = " );
fprintf ( residue_sor, "%d", iter_sor );
fprintf ( residue_sor, "\n" );
/* Data */
for ( i = 1; i <= iter_sor; i++ )
	{
		fprintf (residue_sor, "%d", i);
		fprintf (residue_sor, " " );
		fprintf (residue_sor, "%e", log10(L_infty[i]) );
		fprintf (residue_sor, "\n");
	}
fclose(residue_sor);

/* Error Norms Outputs */
norms_sor=fopen("norms_sor.txt","w");
fprintf ( norms_sor, "Number of Iteration for the Point Gauss-Seidle SOR method = " );
fprintf ( norms_sor,"%d", iter_sor );
fprintf ( norms_sor, "\n" );
fprintf ( norms_sor, "CPU Time (ms) = " );
fprintf ( norms_sor, "%lf", cpu_time_sor );
fprintf ( norms_sor, "\n" );
fprintf ( norms_sor, "L_infty = " );
fprintf ( norms_sor, "%e", L_inf );
fprintf ( norms_sor, "\n");
fprintf ( norms_sor, "L_1 = " );
fprintf ( norms_sor, "%e", L1 );
fprintf ( norms_sor, "\n");
fprintf ( norms_sor, "L_2 = " );
fprintf ( norms_sor, "%e", L2 );
fclose(norms_sor);


printf("Numerical solution using Point Gauss-Seidle SOR with Relaxation Parameter (Omega) of %f has finished.",omega);
printf(" Number of iterations for this method: %d\n\n",iter_sor);


}


void main ()
{

FILE *dat_anal;



// Calculating Space Discretization Parameters
dr = (r_out - r_in)/(imax-1 );
dth = (th_end - th_start)/(jmax-1 );

int i,j;
for ( i = 0; i < imax ; i++)
{
	for ( j = 0; j < jmax ; j++)
	{
		r  [ i ] = r_in     + i * dr ;
		th [ j ] = th_start + j * dth;
		
		x [ i ] [ j ] = r [ i ] * cos ( th [ j ] );
		y [ i ] [ j ] = r [ i ] * sin ( th [ j ] );
	}
}
///////////////////////////////////////////////////////////

// Analytical Solution////////////////////////////////////

////The integ[n] array stores the results of that nasty little definite integral in the analytical solution for n=1 to n=20.

double integ[20] = { -0.33882698, 0.04924979, -0.00518703, 0.00261378, -0.00052642, 0.00059195, -0.00014258, 0.00022581, -0.00005844, 0.00011025,
						-0.00002973, 0.00006217, -0.00001724, 0.00003854, -0.00001091, 0.00002556, -0.00000735, 0.00001782, -0.00000520, 0.00001293 };

double A, sigma, n;

for ( i = 0 ; i < imax ; i++ )
{
	for ( j = 0 ; j < jmax ; j++ )
	{
		sigma = 0.0;

		for ( int n = 1 ; n <=20 ; n++)
		{

			A = - 2*sin(n*M_PI*log(r[i])/log(2))*exp(n*pow(M_PI,2)/(2*log(2)))*integ[n-1]*(exp(n*M_PI*th[j]/log(2)) - exp(-n*M_PI*th[j]/log(2)))/(log(2)*(exp(n*pow(M_PI,2)/log(2)) - 1));
			sigma = sigma + A;
		}

		T_anal [i][j] = (4*sin(2*th[j])*(pow(r[i],4) - 1))/(15*pow(r[i],2)) + sigma;

	}

}


/* Tecplot Output for The Analytical Solution*/
/* Header */
dat_anal=fopen("dat_anal.dat","w");
fprintf ( dat_anal,"Variables = " );
fprintf ( dat_anal, "x");
fprintf ( dat_anal, "," );
fprintf ( dat_anal, "y");
fprintf ( dat_anal, "," );
fprintf ( dat_anal, "T");
fprintf ( dat_anal, "\n");
fprintf ( dat_anal, "Zone I = " );
fprintf ( dat_anal, "%d", imax );
fprintf ( dat_anal, ", J = " );
fprintf ( dat_anal, "%d", jmax );
fprintf ( dat_anal, "\n" );
/* Data */
for ( i = 0; i < imax; i++ )
{
	for ( j = 0; j < jmax; j++ )
	{
		fprintf (dat_anal, "%e", x [i][j]);
		fprintf (dat_anal, " " );
		fprintf (dat_anal, "%e", y [i][j]);
		fprintf (dat_anal, " " );
		fprintf (dat_anal, "%e", T_anal [i][j]);
		fprintf (dat_anal, "\n");
	}
}
fclose(dat_anal);

initialize_all_solutions ();

run_jacubi_method ();
run_gauss_method ();
run_rsweep_method ();
run_thsweep_method ();
run_adi_rth_method ();
run_adi_thr_method ();
run_sor_method ();


printf("Press Enter to terminate the program.\n");
cin.get();

}
