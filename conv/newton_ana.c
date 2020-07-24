
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct {

    int M;				///     'order' of scheme
    int N;				/// 	number of cells on interval
    double xRange;		/// 	xRange of interval
    double xMin;		/// 	smallest x value
    double dx;			///     cell size
    int NSC;            // #subcells

    double CFL;			///     CFL factor

    double tEND;		/// 	final time

// iteration scheme variables - for all problems identical
    int iterMAX;		///     maximum number of iterations --> M+1
    double eps;			/// 	exit condition for iteration scheme

    char project[255];  /// 	Problem definition
    char shape[255];    /// 	initial values
    char scheme[255];   ///     scheme

// convergence test
    int conv_N_start;   ///     lowest resolution
    int conv_N_end;     ///     highest resolution
    int conv_M_start;   ///     lowest order
    int conv_M_end;     ///     highest order
    int change_N;       ///     change factor
    int change_M;       ///     change summand

} PARA;



void Read_Parameter ( char *name, PARA *par );
void Setup_Convergence_Data_Structure ( PARA *par );
void Calc_Analytic_Solution ( double *P, PARA *par );
void Calc_Cell_Average ( double *P, double *R, PARA *par );
void Output_Analytic_Solution ( double *R, PARA *par );
double QiuShu ( double x );
double F_U ( double u, double x, double t );
double F_strich_U ( double u, double x, double t );

// Gauss Legendre nodes on unit interval ...
double LAM[10] = { 0.0130467357414141399610179939578,
                        0.0674683166555077446339516557883,
                        0.160295215850487796882836317443,
                        0.283302302935376404600367028417,
                        0.425562830509184394557586999435,
                        0.574437169490815605442413000565,
                        0.716697697064623595399632971583,
                        0.839704784149512203117163682557,
                        0.932531683344492255366048344212,
                        0.986953264258585860038982006042 };

//... with corresponding weights
double W[10] = {   0.0666713443086881,
                    0.1494513491505806,
                    0.2190863625159820,
                    0.2692667193099963,
                    0.2955242247147529,
                    0.2955242247147529,
                    0.2692667193099963,
                    0.2190863625159820,
                    0.1494513491505806,
                    0.0666713443086881 };




/// ########################################################
/// ###### MAIN ############################################
/// ########################################################


int main ( int argc, char *argv[] ) {

    PARA par;
    Read_Parameter ( argv[1], &par );

    for ( par.M = par.conv_M_start; par.M <= par.conv_M_end; par.M += par.change_M ) {
        par.NSC = 2*par.M + 1;

        for ( par.N = par.NSC*par.conv_N_start; par.N <= par.NSC*par.conv_N_end; par.N *= par.change_N ) {

            Setup_Convergence_Data_Structure ( &par );
            par.dx = (double)(par.xRange) / (double)(par.N);

            double *Points = calloc ( 10*par.N, sizeof ( double ) );
            double *Result = calloc ( par.N, sizeof ( double ) );

            Calc_Analytic_Solution ( Points, &par );
            Calc_Cell_Average ( Points, Result, &par );
            Output_Analytic_Solution ( Result, &par );

            printf ( "\n\n" );

            free ( Points );
            free ( Result );
        }
    }

    return 0;
}


/// ########################################################
/// ########################################################
/// ########################################################



// ### ### ### ### ### ### ###

void Read_Parameter ( char *name, PARA *par ) {
/*
 * name - filename of parameter file
 * par  - parameter
 *
 * read parameter file and store parameter in par
 * set dependent parameters
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
        else if ( strcmp ( var, "scheme" ) == 0 ) {
            sprintf ( par->scheme, "%s", wert );
        }
        else if ( strcmp ( var, "N" ) == 0 ) {
            par->N = atoi ( wert );
        }
        else if ( strcmp ( var, "xRange" ) == 0 ) {
            par->xRange = atof ( wert );
        }
        else if ( strcmp ( var, "xMin" ) == 0 ) {
            par->xMin = atof ( wert );
        }
        else if ( strcmp ( var, "CFL_sub" ) == 0 ) {
            par->CFL = atof ( wert );
        }
        else if ( strcmp ( var, "tEND" ) == 0 ) {
            par->tEND = atof ( wert );
        }
        else if ( strcmp ( var, "eps" ) == 0 ) {
            par->eps = atof ( wert );
        }
        else if ( strcmp ( var, "shape" ) == 0 ) {
            sprintf ( par->shape, "%s", wert );
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
        else if ( strcmp ( var, "iterMAX" ) == 0 ) {
            par->iterMAX = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_M_end" ) == 0 ) {
            par->conv_M_end = atoi ( wert );
        }
        else if ( strcmp ( var, "change_M" ) == 0 ) {
            par->change_M = atoi ( wert );
        }
        else if ( strcmp ( var, "change_N" ) == 0 ) {
            par->change_N = atoi ( wert );
        }
        else {
            printf ( "*** Can't read parameter %s\n", var );
        }
    }
    fclose ( fp );

    par->iterMAX = 1000;
    par->dx = (double)(par->xRange) / (double)(par->N);

}

// ### ### ### ### ### ### ###

void Setup_Convergence_Data_Structure ( PARA *par ) {
/*
 * setting up data structure to store results
 */

    printf ( "\nSetting up data structure for convergence:\n" );

    // 'project'_convergence_M'M'
    char command[255];
    sprintf ( command, "mkdir ../data/%s_convergence_M%d", par->project, par->M );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d\n", par->project, par->M );

    // 'project'_convergence_M'M'/'shape'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s", par->project, par->M, par->shape );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s\n", par->project, par->M, par->shape );

    // 'project'_convergence_M'M'/'shape'/'scheme'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s", par->project, par->M, par->shape, par->scheme );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s\n", par->project, par->M, par->shape, par->scheme );

    // 'project'_convergence_M'M'/'shape'/'scheme'/CFL_'CFL'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s/CFL_%.8f", par->project, par->M, par->shape, par->scheme, par->CFL );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f\n", par->project, par->M, par->shape, par->scheme, par->CFL );

    // 'project'_convergence_M'M'/'shape'/'scheme'/CFL_'CFL'/N_'N'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d", par->project, par->M, par->shape, par->scheme, par->CFL, par->N );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d\n", par->project, par->M, par->shape, par->scheme, par->CFL, par->N );

}

// ### ### ### ### ### ### ###

double QiuShu ( double x ) {
    return ( (0.5 + sin( M_PI * x )) );
}

double F_U ( double u, double x, double t ) {
    return ( u - (0.5 + sin( M_PI*(x - u*t) )) );
}

double F_strich_U ( double u, double x, double t ) {
    return ( 1 + M_PI*t*cos( M_PI*(x - u*t) ) );
}

// ### ### ### ### ### ### ###

static double initial_condition ( int i, int j, PARA *par ) {
    return QiuShu ( par->xMin + (i+LAM[j])*par->dx );
}

static double newton_raphson ( double init, int i, int j, PARA *par ) {

    double change = 100;
    int iter = 0;
    double Old = init;
    double New;
    double t = par->tEND;
    double x = par->xMin + (i+LAM[j])*par->dx;

    do {
        New = Old - ( F_U( Old, x, t ) )/( F_strich_U( Old, x, t ) );

        change = fabs ( Old - New );

        Old = New;

        iter++;
    } while ( iter < par->iterMAX && change > par->eps );

    return Old;
}


void Calc_Analytic_Solution ( double *P, PARA *par ) {

    int i, j;
    double init;

    for ( i = 0; i < par->N; ++i ) {
        for ( j = 0; j < 10; ++j ) {
            init = initial_condition ( i, j, par );
            P[i*10 + j] = newton_raphson ( init, i, j, par );
        }
    }

}

// ### ### ### ### ### ### ###

void Output_Analytic_Solution ( double *R, PARA *par ) {

    int i;
    char file[255];

    sprintf ( file, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/ana.txt", par->project, par->M, par->shape, par->scheme, par->CFL, par->N );

    FILE *fp = fopen ( file, "w" );

    for ( i = 0; i < par->N-1; ++i )
        fprintf ( fp, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, R[i] );

    fprintf ( fp, "%.15f %.15f", par->xMin+(par->N-0.5)*par->dx, R[(par->N-1)] );

    fclose ( fp );

    printf ( "# SUCCESS # - output convergence file: %s\n", file );

}

// ### ### ### ### ### ### ###

void Calc_Cell_Average ( double *P, double *R, PARA *par ) {

    int i, j;
    double temp;

    for ( i = 0; i < par->N; ++i ) {

        temp = 0.0;

        for ( j = 0; j < 10; ++j ) {
            temp += W[j] * P[i*10 + j];
        }
        R[i] = 0.5*temp;
    }

}






























