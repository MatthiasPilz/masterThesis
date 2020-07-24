#include "stdio.h"
#include "stdlib.h"
#include "math.h"


typedef struct {

    int M;				/// 'order' of scheme

// grid setup for equidistant discretisation
    int NSC;            // number of subcells
    int N;				/// 	number of cells on interval
    double CFL;			///     CFL number
    double dx;                  // HACK passend fuer TOV stern!
    double dt;
    double tEND;		/// 	final time
    int output;         ///     time steps between output
    int Nt;             ///     number of total time steps

    double xMin;
    double xRange;

// project specific variables
    char project[255];  /// 	problem definition
    char shape[255];    /// 	initial values
    char scheme[255];   ///     scheme
    char boundary[255]; ///     for TOV only...

    int dimension;

// convergence test
    int conv_N_start;   ///     lowest resolution
    int conv_N_end;     ///     highest resolution
    int conv_M_start;   ///     lowest order
    int conv_M_end;     ///     highest order
    int change_N;       ///     change factor
    int change_M;       ///     change summand

} PARA;

void Read_Parameter ( char *name, PARA *par );




int main ( int argc, char *argv[] ) {

    PARA par;
    Read_Parameter( argv[1], &par );

    double norm_L1 = 0.0;
    double norm_L2 = 0.0;
    double norm_Linf = 0.0;
    char init[1024];
    char final[1024];
    char out[1024];
    char diff[1024];
    double temp_x;
    double temp_val;
    int i;
    int Nt_evolve;

    double Log_N;
    double Log_L1;
    double temp;


    // loop orders of the scheme
    for ( par.M = par.conv_M_start; par.M <= par.conv_M_end; par.M += par.change_M ) {
        printf ( "\n\t M = %d\n", par.M );

        par.NSC = 2*par.M + 1;

        // loop different resolutions
        for ( par.N = par.NSC*par.conv_N_start; par.N <= par.NSC*par.conv_N_end; par.N *= par.change_N ) {

            norm_L1 = 0.0;
            norm_L2 = 0.0;
            norm_Linf = 0.0;

            par.dx = par.xRange / par.N;
            par.dt = par.CFL * par.dx;
            par.Nt = (int) ( par.tEND / (par.CFL * par.dx) );



          for ( Nt_evolve = 0; Nt_evolve < par.Nt; Nt_evolve += par.output ) {
            sprintf ( init, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_ana.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N );
            sprintf ( final, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_Nt_%08d.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N, Nt_evolve );

            FILE *fp_ini = fopen ( init, "r" );
            FILE *fp_fin = fopen ( final, "r" );

            double *I = calloc ( par.N, sizeof ( double ) );
            double *F = calloc ( par.N, sizeof ( double ) );
            double *X = calloc ( par.N, sizeof ( double ) );

            i = 0;
            while ( !feof( fp_ini ) ) {		// go through file line by line
                fscanf ( fp_ini, "%lf, %lf", &temp_x, &temp_val);
                /// HACK HIER
                I[i] = temp_val;			// write result in array
                X[i] = temp_x;
                i++;
            }
            fclose ( fp_ini );

            i = 0;
            while ( !feof( fp_fin ) ) {		// go through file line by line
                fscanf ( fp_fin, "%lf, %lf", &temp_x, &temp_val);
                if ( X[i] != temp_x ) printf ( "ERROR === shift in x-values...i = %d\t %lf != %lf\n", i, X[i], temp_x );
                F[i++] = temp_val;			// write result in array
            }
// 		printf ( "final i = %d\n\n", i );
            fclose ( fp_fin );				// close file

            sprintf ( out, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/error_evolution_M%d_N%d.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N, par.M, (int)(par.N/par.NSC) );
            FILE *f_out = fopen ( out, "a" );

            fprintf ( f_out, "\"Time=%lf", par.dt*Nt_evolve );

            for ( i = 0; i < par.N; ++i ) { 				/// L1-norm
                norm_L1 = ( I[i] - F[i] );
                fprintf ( f_out, "%.15f, %.15f\n", par.xMin+(i+0.5)*par.dx, norm_L1 );
            }

            fclose ( f_out );


//            if ( strcmp( par.boundary, "symmetric" ) == 0 ) {
//            /// DIFFERENCE BETWEEN LEFT AND RIGHT SIDE
//                sprintf ( diff, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/error_diff_M%d_N%d.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N, par.M, (int)(par.N/par.NSC) );
//                FILE *f_diff = fopen ( diff, "a" );
//                fprintf ( f_diff, "\"Time=%lf", par.dt*Nt_evolve );
//                for ( i = 0; i < par.N/2-1; ++i ) {
//                    temp = F[i] - F[par.N-1 - i];
//                    fprintf ( f_diff, "%.15f, %.15f\n", par.xMin+(i+0.5)*par.dx, temp );
//                }
//                fclose ( f_diff );
//            }

            free ( I );
            free ( F );
            free ( X );
          }
        }
    }
    printf ( "\a\n" );
    return 0;
}







void Read_Parameter ( char *name, PARA *par ) {
/*
 * name - filename of parameter file
 * par  - parameter
 *
 * read parameter file and store parameter in par (copied from main program)
 */

    char input[256];
    char var[64];
    char wert[64];
    char filename[255];

    sprintf ( filename, "../par/%s", name );

    FILE *fp = fopen ( filename, "r" );

// set parameters from file
    while ( fgets ( input, 256, fp ) ) {
        sscanf ( input, "%s\t=\t%s", var, wert );
        if ( strcmp ( var, "M" ) == 0 ) {
            par->M = atoi ( wert );
        }
        else if ( strcmp ( var, "project" ) == 0 ) {
            sprintf ( par->project, "%s", wert );
        }
        else if ( strcmp ( var, "N" ) == 0 ) {
            par->N = atoi ( wert );
        }
        else if ( strcmp ( var, "CFL_sub" ) == 0 ) {
            par->CFL = atof ( wert );
        }
        else if ( strcmp ( var, "tEND" ) == 0 ) {
            par->tEND = atof ( wert );
        }
        else if ( strcmp ( var, "xMin" ) == 0 ) {
            par->xMin = atof ( wert );
        }
        else if ( strcmp ( var, "xRange" ) == 0 ) {
            par->xRange = atof ( wert );
        }
        else if ( strcmp ( var, "shape" ) == 0 ) {
            sprintf ( par->shape, "%s", wert );
        }
        else if ( strcmp ( var, "dimension" ) == 0 ) {
            par->dimension = atoi ( wert );
        }
        else if ( strcmp ( var, "scheme" ) == 0 ) {
            sprintf ( par->scheme, "%s", wert );
        }
        else if ( strcmp ( var, "boundary" ) == 0 ) {
            sprintf ( par->boundary, "%s", wert );
        }
        else if ( strcmp ( var, "output" ) == 0 ) {
            par->output = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_N_start" ) == 0 ) {
            par->conv_N_start = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_N_end" ) == 0 ) {
            par->conv_N_end = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_M_start" ) == 0 ) {
            par->conv_M_start = atoi ( wert );
        }
        else if ( strcmp ( var, "change_M" ) == 0 ) {
            par->change_M = atoi ( wert );
        }
        else if ( strcmp ( var, "change_N" ) == 0 ) {
            par->change_N = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_M_end" ) == 0 ) {
            par->conv_M_end = atoi ( wert );
        }
        else {
            printf ( "*** Can't read parameter %s\n", var );
        }
    }
    fclose ( fp );

    // HINT: analytic solution is obtained in external program!
    if ( strcmp ( par->project, "burger" ) == 0 )
        printf ( "*** ATTENTION - BURGER's EQUATION!\n*** Make sure to calc analytic solution beforehand!\n" );

    par->NSC = 2*par->M + 1;

        /// distinguish between 1D and 3D TOV evolution
    if ( strcmp ( par->shape, "TOV" ) == 0 )
        sprintf ( par->shape, "%s_%dD", par->shape, par->dimension );

}


