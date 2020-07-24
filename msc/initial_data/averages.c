#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "TOV.h"

// Collection of all possible initial conditions for ADER-DG program


double Ue0 ( double x );
double Ve0 ( double x );
double sgauss ( double x );
double gauss ( double x );
double Dgauss ( double x );
double gaussadd ( double x );
double step ( double x );

double paper ( double x );
double qiushu ( double x );
double monte ( double x );
double monteexp ( double x );

double high_gauss ( double x );
double high_Dgauss ( double x );

double bugner ( double x );
double Dbugner ( double x );

double bugner2 ( double x );
double Dbugner2 ( double x );

double Ibarrier ( double x );
double barrier ( double x );

double Itriangle ( double x );
double triangle ( double x );

double SOD_r ( double x );
double SOD_p ( double x );
double SOD_v ( double x );

double paper_bulev ( double x );

double QS_p ( double x );
double QS_r ( double x );
double QS_v ( double x );

double LAX_p ( double x );
double LAX_r ( double x );
double LAX_v ( double x );

double MACH_p ( double x );
double MACH_r ( double x );
double MACH_v ( double x );

double Blast_p ( double x );
double Blast_r ( double x );
double Blast_v ( double x );

double grhdc_sphere_p ( double x );
double grhdc_sphere_r ( double x );
double grhdc_sphere_v ( double x );
double grhdc_sphere_eps ( double x );
double grhdc_sphere_D ( double x );
double grhdc_sphere_S ( double x );
double grhdc_sphere_tau ( double x );
double grhdc_sphere_psi4 ( double x );
double grhdc_sphere_alpha ( double x );

double grhdc_sinwave_p ( double x );
double grhdc_sinwave_r ( double x );
double grhdc_sinwave_v ( double x );
double grhdc_sinwave_eps ( double x );
double grhdc_sinwave_D ( double x );
double grhdc_sinwave_S ( double x );
double grhdc_sinwave_tau ( double x );
double grhdc_sinwave_psi4 ( double x );
double grhdc_sinwave_alpha ( double x );

int main ( int argc, char *argv[] ) {
    int N;	//number of cells
	double delx;	//uniform grid
    double xmin;
    double xRange;

    double PI = 3.14159265358979323846264338327950288419716939937510582097494459;

	int i,j;
	double a;
	double b;
		
	char filename[255];		// for output
	char project[255];		// equation to solve --> defines what kind of data is necessary
	char shape[255];		// shape of initial profile
	char command[255];		// for setting up the data structure

	
    // Gauss Integration nodes and weights for 11th order exact cell averages
    double LAMBDA[5] = {-0.90617984593866399279762687829939296512565191076253,
                        -0.53846931010568309103631442070020880496728660690556,
                        0.0,
                        0.53846931010568309103631442070020880496728660690556,
                        0.90617984593866399279762687829939296512565191076253 };
							
    double W[5] = {		0.23692688505618908751426404071991736264326000221241,
                        0.47862867049936646804129151483563819291229555334314,
                        0.56888888888888888888888888888888888888888888888889,
                        0.47862867049936646804129151483563819291229555334314,
                        0.23692688505618908751426404071991736264326000221241 };

    int order = 5;
	
	// set important parameter from input
	N = atoi ( argv[1] );
    xmin = atof ( argv[2] );
    xRange = atof ( argv[3] );
    sprintf ( project, "%s", argv[4] );
    sprintf ( shape, "%s", argv[5] );

    double SinStretch = 2.0 * PI / xRange;


	
    printf ( "Create initial data for %d points in [%.3f, %.3f]\n", N, xmin, xmin+xRange );

	// setup data structure
	printf ( "\nSetting up data structure for initial data:\n" );
	
	sprintf ( command, "mkdir ../math/initial_data/%s", project );
	if ( !system ( command ) ) 
		printf ( "\t created path: ../math/initial_data/%s\n", project );
		
	sprintf ( command, "mkdir ../math/initial_data/%s/%s", project, shape );
	if ( !system ( command ) )
		printf ( "\t created path: ../math/initial_data/%s/%s\n", project, shape );
		
	printf ( "\n" );	
	
    delx = (double) xRange / N;

    // if-else structure to differ between different HCL
	if ( strcmp ( project, "wave" ) == 0 ) { 
    // initial data for wave equation needs cell averages and first spatial derivative
		double *result;		// contains cell averages
		result = (double*) calloc ( N, sizeof (double ) );
			
		double *resultD;	// contains spacial derivative
		resultD = (double*) calloc ( N, sizeof ( double ) );
		
		
		if ( strcmp ( shape, "gauss" ) == 0 ) {	// gauss
			for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));
		
				result[i] = 0.0;
                for ( j = 0; j < order; j++ )
					result[i] += ( W[j] * gauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
			
				result[i] *= 0.5;
			}
	
			for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));
		
				resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
					resultD[i] += ( W[j] * Dgauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
		
				resultD[i] *= 0.5;
			}
        }
        else if ( strcmp ( shape, "high_gauss" ) == 0 ) {	// gauss with higher amplitude
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * high_gauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * high_Dgauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "bugner" ) == 0 ) {	// test function by marcus
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * bugner( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * Dbugner( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "bugner2" ) == 0 ) {	// 2nd test function by marcus
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * bugner2( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * Dbugner2( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "barrier" ) == 0 ) {	// barrier in phi comp --> to avoid inf
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * Ibarrier( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * barrier( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "triangle" ) == 0 ) {	// triangle in phi comp --> to avoid inf
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * Itriangle( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * triangle( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "step") == 0 ) { // unnecessary and derivative wrong!
            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));
                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * step( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }
            for ( i = 0; i < N; i++ ) {	//derivative
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));
                resultD[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    resultD[i] += ( W[j] * step( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                resultD[i] *= 0.5;
            }
            printf ( "*** ERROR - derivative of step function not well defined \nBE CAREFUL!\n" );
        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );



    //output to file
        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N; i++ ) {
            fprintf ( fp, "%.18f %.18f\n", delx*(i+0.5), result[i] );
        }
        fclose ( fp );
        free ( result );

        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICD.txt", project, shape, N );
        FILE * fpD;
        fpD = fopen ( filename, "wx" );
        for ( i = 0; i < N; i++ ) {
            fprintf ( fpD, "%.18f %.18f\n", delx*(i+0.5), resultD[i] );
        }
        fclose ( fpD );
        free ( resultD );

    }

    else if ( strcmp ( project, "advection" ) == 0 ) {
        double *result;		// contains cell averages
        result = (double*) calloc ( N, sizeof (double ) );

        if ( strcmp ( shape, "gauss" ) == 0 ) {	// gauss
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * gauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "triangle" ) == 0 ) {	// triangle
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * Itriangle( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "sinus" ) == 0 ) {	// sinus
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * sin( SinStretch*(0.5*((b-a)*LAMBDA[j] + (b+a))) ) );

                result[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "gaussadd" ) == 0 ) {	// gauss + 1.0
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * gaussadd( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }
        }
        else if ( strcmp ( shape, "barrier") == 0 ) {
            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));
                // HACK
                result[i] = barrier ( 0.5*(a+b) );
            }
        }
        else if ( strcmp ( shape, "paper" ) == 0 ) {	// initial data from paper
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * paper( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );


    //output to file
        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N; i++ ) {
            fprintf ( fp, "%.18f %.18f\n", delx*(i+0.5), result[i] );
        }
        fclose ( fp );
        free ( result );

    }

    else if ( strcmp ( project, "burger" ) == 0 ) {
        double *result;		// contains cell averages
        result = (double*) calloc ( N, sizeof (double ) );

        if ( strcmp ( shape, "gauss" ) == 0 ) {	// gauss
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * gauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "sgauss" ) == 0 ) {	// small gauss
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * sgauss( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "qiushu" ) == 0 ) {	// paper qiu shu
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * qiushu( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "montelin" ) == 0 ) {	// initial data from paper
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * paper( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "montesquare" ) == 0 ) {	// initial data from paper
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * monte( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else if ( strcmp ( shape, "monteexp" ) == 0 ) {	// initial data from paper
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * monteexp( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );


    //output to file
        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N; i++ ) {
            fprintf ( fp, "%.18f %.18f\n", delx*(i+0.5), result[i] );
        }
        fclose ( fp );
        free ( result );

    }

    else if ( strcmp ( project, "grhdc" ) == 0 ) {
        double *result;		// contains cell averages
        result = (double*) calloc ( N*9, sizeof (double ) );

        if ( strcmp ( shape, "sphere" ) == 0 ) {	// gauss
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ ) {
                    result[i + 0*N] += 0.5*( W[j] * grhdc_sphere_D    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //D
                    result[i + 1*N] += 0.5*( W[j] * grhdc_sphere_S    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //S
                    result[i + 2*N] += 0.5*( W[j] * grhdc_sphere_tau  ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //TAU

                    result[i + 3*N] += 0.5*( W[j] * grhdc_sphere_alpha( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //ALPHA
                    result[i + 4*N] += 0.5*( W[j] * grhdc_sphere_psi4 ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //PSI4

                    result[i + 5*N] += 0.5*( W[j] * grhdc_sphere_r    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //RHO
                    result[i + 6*N] += 0.5*( W[j] * grhdc_sphere_v    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //V
                    result[i + 7*N] += 0.5*( W[j] * grhdc_sphere_eps  ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //EPS
                    result[i + 8*N] += 0.5*( W[j] * grhdc_sphere_p    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //P
                }

            }

        }

        else if ( strcmp ( shape, "sinwave" ) == 0 ) {	
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ ) {
                    result[i + 0*N] += 0.5*( W[j] * grhdc_sinwave_D    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //D
                    result[i + 1*N] += 0.5*( W[j] * grhdc_sinwave_S    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //S
                    result[i + 2*N] += 0.5*( W[j] * grhdc_sinwave_tau  ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //TAU

                    result[i + 3*N] += 0.5*( W[j] * grhdc_sinwave_alpha( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //ALPHA
                    result[i + 4*N] += 0.5*( W[j] * grhdc_sinwave_psi4 ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //PSI4

                    result[i + 5*N] += 0.5*( W[j] * grhdc_sinwave_r    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //RHO
                    result[i + 6*N] += 0.5*( W[j] * grhdc_sinwave_v    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //V
                    result[i + 7*N] += 0.5*( W[j] * grhdc_sinwave_eps  ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //EPS
                    result[i + 8*N] += 0.5*( W[j] * grhdc_sinwave_p    ( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) ); //P
                }

            }

        }
        
        else if ( (strcmp ( shape, "TOV_3D" ) == 0) || (strcmp ( shape, "TOV_1D" ) == 0) ) {

        	int i,j;
			static double *p_r,*p_m,*p_h,*p_rho,*p_pre,*p_phi,*p_riso;
			int npts = 10000;

			double rhoc = 0.00128;
			double R0   = 10.;
			double perturb  = 0.;
            double atm    = 1e-8;
			double atmfac = 100.;

            double Gamma  = 2.;
            double kappa  = 100.;

			double rhoatm  = atm * rhoc;
			double rhoatmf = rhoatm * atmfac;
			double patm    = kappa * pow(rhoatm,Gamma);
			double eatm    = kappa/(Gamma-1.)*pow(rhoatm,Gamma-1.);
		
			/* size of the array is not yet known -> is handled automatically 
			 prototypes are different ! */

			tov_r(rhoc, R0, &npts, &p_r, &p_m, &p_rho, 
				   &p_pre, &p_phi, &p_riso);

			double M = p_m[npts-1];
			double R = p_riso[npts-1];
			double r,phi,dphi,rsch,drsch,e,pres,dummy,p,dpdrho,rho;    
            double grho, gp, geps, gD, gS, gtau, gpsi4;    		   

			// set TOV to grid
            for ( i = 0; i < N; i++ ) {	//cell averages
              a = xmin + (delx * i);
              b = xmin + (delx * (i+1));

              for ( j = 0; j < order; j++ ) {			  
                 
              r = 0.5*((b-a)*LAMBDA[j] + (b+a));
			  r = sqrt(r*r);

			  if (r<R) {
				interp_lag4(p_r, p_riso, npts, r, &rsch, &drsch,  &dummy);
				interp_lag4(p_phi, p_riso, npts, r, &phi,&dphi,  &dummy);
                gpsi4 = pow(rsch/r,2.);
				result[i + 3*N] += 0.5*( W[j] * exp(phi) ); //ALPHA
				result[i + 4*N] += 0.5*( W[j] * gpsi4 );//PSI4
				interp_lag4(p_rho, p_riso, npts, r, &grho,  &dummy, &dummy);
				eos_simple(grho, &gp, &dpdrho, &geps);
				if ( r==0 ) {
				r=1.e-13;
				interp_lag4(p_r, p_riso, npts, r, &rsch, &drsch,  &dummy);
				result[i + 4*N] += 0.5*( W[j] * pow(rsch/r,2.) ); //PSI4
			    }
                result[i + 5*N] += 0.5*( W[j] * grho ); //RHO
                result[i + 6*N] += 0.5*( W[j] * 0.0  ); //V
                result[i + 7*N] += 0.5*( W[j] * geps ); //EPS
                result[i + 8*N] += 0.5*( W[j] * gp   ); //P
			  } else {
                grho = 0.0;
                gpsi4 = pow(1.0+0.5*M/r, 4.);
                eos_simple(grho, &gp, &dpdrho, &geps);
				result[i + 3*N] += 0.5*( W[j] * (1.0-0.5*M/r)/(1.0+0.5*M/r) );//ALPHA
				result[i + 4*N] += 0.5*( W[j] * gpsi4); //PSI4
                result[i + 5*N] += 0.5*( W[j] * grho ); //RHO
                result[i + 7*N] += 0.5*( W[j] * geps ); //EPS
                result[i + 8*N] += 0.5*( W[j] * gp   ); //P
			  }
			  			  			  
			 /* set atm*/
			 if (grho<rhoatmf){
               grho = rhoatm;
               geps = eatm;
			   result[i + 5*N] += 0.5*( W[j] * rhoatm );
			   result[i + 7*N] += 0.5*( W[j] * eatm );
			 }

			 // conservative variables
			 grhdc_p2c_point(grho, geps, 0.0, 0.0, 0.0, gp, gpsi4,
			                &gD, &gtau, &gS, &dummy,  &dummy);		
			
             result[i + 0*N] += 0.5*( W[j] * gD  ); //D
             result[i + 1*N] += 0.5*( W[j] * gS  ); //S
             result[i + 2*N] += 0.5*( W[j] * gtau); //TAU

	        } //end loop over gridpoints 
          }

        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );


    //output to file
        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N; i++ ) {
            fprintf ( fp, "%.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f\n", delx*(i+0.5), 
                      result[i + 0*N], result[i + 1*N], result[i + 2*N], result[i + 3*N], 
                      result[i + 4*N], result[i + 5*N], result[i + 6*N], result[i + 7*N], result[i + 8*N] );
        }
        fclose ( fp );
        free ( result );

    }

    else if ( strcmp ( project, "bulev" ) == 0 ) {
        double *result;		// contains cell averages
        result = (double*) calloc ( N, sizeof (double ) );

        if ( strcmp ( shape, "paper" ) == 0 ) {	// initial data from paper
            for ( i = 0; i < N; i++ ) {	//cell averages
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                result[i] = 0.0;
                for ( j = 0; j < order; j++ )
                    result[i] += ( W[j] * paper_bulev( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                result[i] *= 0.5;
            }

        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );


    //output to file
        sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N; i++ ) {
            fprintf ( fp, "%.18f %.18f\n", delx*(i+0.5), result[i] );
        }
        fclose ( fp );
        free ( result );

    }

    else if ( strcmp ( project, "euler" ) == 0 ) {
        double *pressure;
        pressure = calloc ( N, sizeof ( double ) );

        double *density;
        density = calloc ( N, sizeof ( double ) );

        double *velocity;
        velocity = calloc ( N, sizeof ( double ) );

        if ( strcmp ( shape, "SOD" ) == 0 ) {
            printf ( "*** ATTENTION! - discontinuity at x=0.5\n\tadjust xmin and xRange accordingly!\n\n" );

            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                pressure[i] = 0.0;
                density[i]  = 0.0;
                velocity[i] = 0.0;

                for ( j = 0; j < order; j++ ) {
                    pressure[i] += ( W[j] * SOD_p( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    density[i]  += ( W[j] * SOD_r( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    velocity[i] += ( W[j] * SOD_v( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                }

                pressure[i] *= 0.5;
                density[i]  *= 0.5;
                velocity[i] *= 0.5;

            }

        }
        else if ( strcmp ( shape, "qiushu" ) == 0 ) {

            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                density[i]  = 0.0;

                for ( j = 0; j < order; j++ )
                    density[i]  += ( W[j] * QS_r( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                pressure[i] = 1.0;
                density[i]  *= 0.5;
                velocity[i] = 0.7;

            }

        }
        else if ( strcmp ( shape, "LAX" ) == 0 ) {
            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                pressure[i] = 0.0;
                density[i]  = 0.0;
                velocity[i] = 0.0;

                for ( j = 0; j < order; j++ ) {
                    pressure[i] += ( W[j] * LAX_p( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    density[i]  += ( W[j] * LAX_r( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    velocity[i] += ( W[j] * LAX_v( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                }

                pressure[i] *= 0.5;
                density[i]  *= 0.5;
                velocity[i] *= 0.5;

            }
        }
        else if ( strcmp ( shape, "MACH" ) == 0 ) {
            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                pressure[i] = 0.0;
                density[i]  = 0.0;
                velocity[i] = 0.0;

                for ( j = 0; j < order; j++ ) {
                    pressure[i] += ( W[j] * MACH_p( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    density[i]  += ( W[j] * MACH_r( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );
                    velocity[i] += ( W[j] * MACH_v( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                }

                pressure[i] *= 0.5;
                density[i]  *= 0.5;
                velocity[i] *= 0.5;

            }
        }
        else if ( strcmp ( shape, "Blast" ) == 0 ) {
            for ( i = 0; i < N; ++i ) {
                a = xmin + (delx * i);
                b = xmin + (delx * (i+1));

                pressure[i] = 0.0;

                for ( j = 0; j < order; j++ ) {
                    pressure[i] += ( W[j] * Blast_p( 0.5*((b-a)*LAMBDA[j] + (b+a)) ) );

                }

                pressure[i] *= 0.5;
                density[i]  = 1.0;
                velocity[i] = 0.0;

            }
        }
        else
            printf ( "*** ERROR - initial data - unknown shape!\n" );

       //output to file
       /// pressure
        sprintf ( filename, "../math/initial_data/%s/%s/%d_pressure.txt", project, shape, N );
        FILE * fp;
        fp = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N-1; i++ ) {
            fprintf ( fp, "%.18f %.18f\n", delx*(i+0.5), pressure[i] );
        }
        fprintf ( fp, "%.18f %.18f\n", delx*(N-0.5), pressure[N-1] );
        fclose ( fp );

        /// density
        sprintf ( filename, "../math/initial_data/%s/%s/%d_density.txt", project, shape, N );
        FILE * fp1;
        fp1 = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N-1; i++ ) {
            fprintf ( fp1, "%.18f %.18f\n", delx*(i+0.5), density[i] );
        }
        fprintf ( fp1, "%.18f %.18f\n", delx*(N-0.5), density[N-1] );
        fclose ( fp1 );

        /// velocity
        sprintf ( filename, "../math/initial_data/%s/%s/%d_velocity.txt", project, shape, N );
        FILE * fp2;
        fp2 = fopen ( filename, "wx" );  // !!! wx -> segfault if file exists!
        for ( i = 0; i < N-1; i++ ) {
            fprintf ( fp2, "%.18f %.18f\n", delx*(i+0.5), velocity[i] );
        }
        fprintf ( fp2, "%.18f %.18f\n", delx*(N-0.5), velocity[N-1] );
        fclose ( fp2 );

    }


    else
        printf ( "*** ERROR - initial data - unknown project!\n" );

	

	



		
	return 0;

}


// #### #### ####

double Ue0 ( double x ) {
	double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
	return ( 4.0 + 0.1*sin( 2*PI*x ) );
}

double Ve0 ( double x ) {
	double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
	return ( 6.0 + 0.3*cos( 2*PI*x ) );
}

// #### #### ####

double paper ( double x ) {
    if ( x <= 0.1 )
        return 2.0;
    else if ( x >= 0.1 && x <= 0.2 )
        return ( 2.0 + (x - 0.1)*10.0 );
    else if ( x >= 0.2 && x <= 0.4 )
        return 3.0;
    else if ( x >= 0.4 && x <= 0.6 )
        return ( 3.0 - (x - 0.4)*5.0 );
    else if ( x >= 0.6 )
        return 2.0;
    else
        printf ( "*** ERROR - x range\n" );

}

double qiushu ( double x ) {
	double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
    return ( 0.5 + sin( PI * x ) );

}

double monte ( double x ) {
    double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
    return ( sin ( 2*PI*x ) );
}

static double omega ( double x ) {
    return ( (0.3 - x)/(sqrt((0.3-x)*(0.3-x) + 1e-4)) );
}

double monteexp ( double x ) {
    return ( (1-omega(x)) + ((1+omega(x))/2.0) );
}

// #### #### ####

double sgauss ( double x ) {
    return ( exp ( -pow((x - 0.5),2) * 100.0 ) );
}

// #### #### ####

double gauss ( double x ) {
	return ( 2.0*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

double Dgauss ( double x ) {
    return ( -600*(x-0.5)*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

double gaussadd ( double x ) {
    return ( 1.0 + 2.0*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

// #### #### ####

double step ( double x ) {
    if ( x < 0.5 )
        return 1.0;
    else
        return 0.0;
}

// #### #### ####

double high_gauss ( double x ) {
    return ( 5.0*exp ( -pow( (x-0.5), 2 ) * 300.0 ) );
}

double high_Dgauss (double x ) {
    return ( -3000*(x-0.5)*exp( -pow( (x-0.5), 2 ) * 300.0 ) );
}

// #### #### ####

double bugner ( double x ) {
    return ( -30.0*(x-0.5)*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

double Dbugner ( double x ) {
    return ( ( 2220.0 - 9000.0*x + 9000.0*x*x )*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

// #### #### ####

double bugner2 ( double x ) {
    return ( 100.0 / pow(1.0 + 1.0/((-1.0 + x)*x),4) );
}

double Dbugner2 ( double x ) {
    return ( 400.0*pow((x-1.0),3)*(x*x*x)*(2.0*x-1.0)/(pow((1-x+x*x),5)) );
}

// #### #### ####

double barrier ( double x ) {

    if ( x < 0.35 )
        return 0.0;
    else if ( x >= 0.35 && x < 0.65 )
        return 1.0;
    else if ( x >= 0.65 )
        return 0.0;
    else {
        printf ( "ERROR-average barrier\tx=%lf\n", x );
        return 42;
    }
}

double Ibarrier ( double x ) {

    if ( x < 0.35 )
        return 0.0;
    else if ( x >= 0.35 && x < 0.65 )
        return (x-0.35);
    else if ( x >= 0.65 )
        return 0.3;
    else {
        printf ( "ERROR-average barrier\tx=%lf\n", x );
        return 42;
    }

}

// #### #### ####

double Itriangle ( double x ) {

    if ( x < 0.3 )
        return 0.0;
    else if ( x >= 0.3 && x < 0.5 )
        return (  10*x - 3 );
    else if ( x >= 0.5 && x < 0.7 )
        return ( -10*x + 7 );
    else if ( x >= 0.7 )
        return 0.0;
    else {
        printf ( "ERROR-average triangle\tx=%lf\n", x );
        return 42;
    }

}

double triangle ( double x ) {

    if ( x < 0.3 )
        return 0.0;
    else if ( x >= 0.3 && x < 0.5 )
        return 1.0;
    else if ( x >= 0.5 && x < 0.7 )
        return -1.0;
    else if ( x >= 0.7 )
        return 0.0;
    else {
        printf ( "ERROR-average triangle\tx=%lf\n", x );
        return 42;
    }

}

// #### #### ####

double SOD_r ( double x ) {
// hardcoded discontinuity at x=0.5

    if ( x <= 0.5 )
        return 1.0;
    else
        return 0.125;

}

double SOD_p ( double x ) {
// hardcoded discontinuity at x=0.5

    if ( x <= 0.5 )
        return 1.0;
    else
        return 0.1;

}

double SOD_v ( double x ) {
// hardcoded discontinuity at x=0.5

    if ( x <= 0.5 )
        return 0.0;
    else
        return 0.0;

}

// #### #### ####

double paper_bulev ( double x ) {

    if ( x <= 0 && x >= -0.5 ) {
        return 1.0;
    }
    else
        return 0.0;

}


// ########################################
// ########################################

double QS_p ( double x ) {
    return 1.0;
}

double QS_r ( double x ) {
    double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
    return ( 1.0 + 0.2*sin( PI * x ) );
}

double QS_v ( double x ) {
    return 0.3;
}

// #### #### ####

double LAX_r ( double x ) {

    if ( x <= 0.0 )
        return 0.445;
    else
        return 0.5;

}

double LAX_p ( double x ) {

    if ( x <= 0.0 )
        return 3.528;
    else
        return 0.571;

}

double LAX_v ( double x ) {

    if ( x <= 0.0 )
        return 0.698;
    else
        return 0.0;

}

// #### #### ####

double MACH_r ( double x ) {

    if ( x <= -4.0 )
        return 3.857143;
    else
        return ( 1.0 + 0.2*sin( 5 * x ));

}

double MACH_p ( double x ) {

    if ( x <= -4.0 )
        return 10.3333333;
    else
        return 1.0;

}

double MACH_v ( double x ) {

    if ( x <= -4.0 )
        return 2.629369;
    else
        return 0.0;

}

// #### #### ####

double Blast_r ( double x ) {

}

double Blast_p ( double x ) {

    if ( x < 0.1 )
        return 10.0;
    else if ( x >= 0.1 && x < 0.9 )
        return 1;
    else if ( x >= 0.9 )
        return 5.0;
    else
        printf ( "ERROR - BLAST\n" );

}

double Blast_v ( double x ) {

}

// ########################################
// ########################################

double grhdc_sphere_p ( double x ) {
    return pow(grhdc_sphere_r ( x ), 5./3.);
}

double grhdc_sphere_r ( double x ) {
    return exp(-x*x/(0.25*0.25)) +1.;
}

double grhdc_sphere_v ( double x ) {
    return 0.0;
}

double grhdc_sphere_eps ( double x ) {
    return 1.5 * grhdc_sphere_p ( x ) / grhdc_sphere_r ( x );
}

double grhdc_sphere_D ( double x ) {
double W     = 1/sqrt(1 - grhdc_sphere_psi4(x)*grhdc_sphere_v(x)*grhdc_sphere_v(x));
double h     = 1 + grhdc_sphere_eps(x) + grhdc_sphere_p(x)/grhdc_sphere_r(x);
double gam   = pow(grhdc_sphere_psi4(x),1.5);
double gWWrh = gam * grhdc_sphere_r(x) *W*W * h;
    return gam * grhdc_sphere_r(x) * W;
}

double grhdc_sphere_S ( double x ) {   
double W     = 1/sqrt(1 - grhdc_sphere_psi4(x)*grhdc_sphere_v(x)*grhdc_sphere_v(x));
double h     = 1 + grhdc_sphere_eps(x) + grhdc_sphere_p(x)/grhdc_sphere_r(x);
double gam   = pow(grhdc_sphere_psi4(x),1.5);
double gWWrh = gam * grhdc_sphere_r(x) *W*W * h; 
    return gWWrh * grhdc_sphere_psi4(x) * grhdc_sphere_v(x);
}

double grhdc_sphere_tau ( double x ) {
double W     = 1/sqrt(1 - grhdc_sphere_psi4(x)*grhdc_sphere_v(x)*grhdc_sphere_v(x));
double h     = 1 + grhdc_sphere_eps(x) + grhdc_sphere_p(x)/grhdc_sphere_r(x);
double gam   = pow(grhdc_sphere_psi4(x),1.5);
double gWWrh = gam * grhdc_sphere_r(x) *W*W * h; 
    return gWWrh - gam * (  grhdc_sphere_p(x) + grhdc_sphere_r(x) * W);
}

double grhdc_sphere_psi4 ( double x ) {
    return 1.0;
}

double grhdc_sphere_alpha ( double x ) {
    return 1.0;
}

/**********************************************************************************************/

double grhdc_sinwave_p ( double x ) {
    return 1.0;
}

double grhdc_sinwave_r ( double x ) {
    double PI = 3.14159265358979323846264338327950288419716939937510582097494459;
    return 1.0 + 0.2 * sin(x*PI);
}

double grhdc_sinwave_v ( double x ) {
    return 0.2;
}

double grhdc_sinwave_eps ( double x ) {
    return grhdc_sinwave_p ( x ) / grhdc_sinwave_r ( x );
}

double grhdc_sinwave_D ( double x ) {
double W     = 1/sqrt(1 - grhdc_sinwave_psi4(x)*grhdc_sinwave_v(x)*grhdc_sinwave_v(x));
double h     = 1 + grhdc_sinwave_eps(x) + grhdc_sinwave_p(x)/grhdc_sinwave_r(x);
double gam   = pow(grhdc_sinwave_psi4(x),1.5);
double gWWrh = gam * grhdc_sinwave_r(x) *W*W * h;
    return gam * grhdc_sinwave_r(x) * W;
}

double grhdc_sinwave_S ( double x ) {   
double W     = 1/sqrt(1 - grhdc_sinwave_psi4(x)*grhdc_sinwave_v(x)*grhdc_sinwave_v(x));
double h     = 1 + grhdc_sinwave_eps(x) + grhdc_sinwave_p(x)/grhdc_sinwave_r(x);
double gam   = pow(grhdc_sinwave_psi4(x),1.5);
double gWWrh = gam * grhdc_sinwave_r(x) *W*W * h; 
    return gWWrh * grhdc_sinwave_psi4(x) * grhdc_sinwave_v(x);
}

double grhdc_sinwave_tau ( double x ) {
double W     = 1/sqrt(1 - grhdc_sinwave_psi4(x)*grhdc_sinwave_v(x)*grhdc_sinwave_v(x));
double h     = 1 + grhdc_sinwave_eps(x) + grhdc_sinwave_p(x)/grhdc_sinwave_r(x);
double gam   = pow(grhdc_sinwave_psi4(x),1.5);
double gWWrh = gam * grhdc_sinwave_r(x) *W*W * h; 
    return gWWrh - gam * (  grhdc_sinwave_p(x) + grhdc_sinwave_r(x) * W);
}

double grhdc_sinwave_psi4 ( double x ) {
    return 1.0;
}

double grhdc_sinwave_alpha ( double x ) {
    return 1.0;
}











