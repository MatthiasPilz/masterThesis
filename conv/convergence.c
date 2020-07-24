#include "stdio.h"
#include "stdlib.h"
#include "math.h"


typedef struct {

    int M;				/// 'order' of scheme

// grid setup for equidistant discretisation
    int NSC;            // number of subcells
    int N;				/// 	number of cells on interval
    double CFL;			///     CFL number
    double tEND;		/// 	final time

// project specific variables
    char project[255];  /// 	problem definition
    char shape[255];    /// 	initial values
    char scheme[255];   ///     scheme

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

    double R1;  // relativ outer radius
    if ( argc > 2 && strcmp(par.project, "grhdc") == 0 )
        R1 = atof ( argv[2] );
    else
        R1 = INFINITY;


    double norm_L1 = 0.0;
    double norm_L2 = 0.0;
    double norm_Linf = 0.0;
    char init[1024];
    char final[1024];
    double temp_x;
    double temp_val;
    int i;
    int j;

    double Log_N;
    double Log_L1;


    // loop orders of the scheme
    for ( par.M = par.conv_M_start; par.M <= par.conv_M_end; par.M += par.change_M ) {
        printf ( "\n\t M = %d\n", par.M );

        par.NSC = 2*par.M + 1;

        // loop different resolutions
        for ( par.N = par.NSC*par.conv_N_start; par.N <= par.NSC*par.conv_N_end; par.N *= par.change_N ) {

            norm_L1 = 0.0;
            norm_L2 = 0.0;
            norm_Linf = 0.0;

            sprintf ( init, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_ana.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N );
            sprintf ( final, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_final.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N );

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


            for ( i = par.M; i < par.N-par.M; ++i ) 				/// L2-norm
                norm_L2 += pow ( I[i] - F[i], 2 );

            norm_L2 = norm_L2 / ( par.N-2*par.M );
            norm_L2 = sqrt ( norm_L2 );

            j = 0;
            for ( i = par.M; i < par.N-par.M; ++i ) 				/// L1-norm
                if ( fabs(X[i]) < R1 ) {
                    norm_L1 += fabs ( I[i] - F[i] );
                    j++;
                }

            norm_L1 = norm_L1 / ( j );

            for ( i = par.M; i < par.N-par.M; ++i )                 /// Linf-norm
                if ( fabs(X[i]) < R1 )
                    if ( fabs( I[i] - F[i] ) > norm_Linf )
                        norm_Linf = fabs ( I[i] - F[i] );

            Log_N = log10 ( par.N );
            Log_L1 = log10 ( norm_L1 );

            printf ( "%d\t %.2e\t%.2e\t\t\n", par.N, norm_L1, norm_Linf );

            free ( I );
            free ( F );
            free ( X );
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
        else if ( strcmp ( var, "shape" ) == 0 ) {
            sprintf ( par->shape, "%s", wert );
        }
        else if ( strcmp ( var, "dimension" ) == 0 ) {
            par->dimension = atoi ( wert );
        }
        else if ( strcmp ( var, "scheme" ) == 0 ) {
            sprintf ( par->scheme, "%s", wert );
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
