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
    double temp_x;
    int Nt_evolve;

    int i;

    double I;
    double F;
    double X;


    // loop orders of the scheme
    for ( par.M = par.conv_M_start; par.M <= par.conv_M_end; par.M += par.change_M ) {
        printf ( "\n\t M = %d\n", par.M );

        par.NSC = 2*par.M + 1;

        // loop different resolutions
        for ( par.N = par.NSC*par.conv_N_start; par.N <= par.NSC*par.conv_N_end; par.N *= par.change_N ) {

            norm_L1 = 0.0;
            norm_L2 = 0.0;
            norm_Linf = 0.0;
            par.dx = 24.0 / par.N;  // HACK!!!
            par.dt = par.CFL * par.dx;
            par.Nt = (int) ( par.tEND / (par.CFL * par.dx) );



            sprintf ( out, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/error_origin_M%d_N%d.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N, par.M, par.N );
            FILE *f_out = fopen ( out, "w" );

            for ( Nt_evolve = 0; Nt_evolve < par.Nt; Nt_evolve += par.output ) {
                sprintf ( init, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_ana.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N );
                sprintf ( final, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_Nt_%08d.txt", par.project, par.M, par.shape, par.scheme, par.CFL, par.N, Nt_evolve );

                FILE *fp_ini = fopen ( init, "r" );
                FILE *fp_fin = fopen ( final, "r" );

                for ( i = 0; i < (int)(par.N/2); ++i )
                    fscanf ( fp_ini, "%lf, %lf", &X, &I );
                fclose ( fp_ini );

                for ( i = 0; i < (int)(par.N/2); ++i )
                    fscanf ( fp_fin, "%lf, %lf", &temp_x, &F );
                if ( X != temp_x ) printf ( "ERROR === shift in x-values... %lf != %lf\n", X, temp_x );
                fclose ( fp_fin );				// close file

                norm_L1 = I - F;
                fprintf ( f_out, "%.15f, %.15f\n", Nt_evolve*par.dt, norm_L1 );

            }

            fclose ( f_out );
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


