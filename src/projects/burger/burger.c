#include "ader.h"
#include "burger.h"

void Burger_Var_Num ( PARA *par ) {

    par->var = 1;
    printf ( "burger equation needs %d variables\n", par->var );

}

// ### ### ### ### ### ### ###

static void create_initial_data_burger ( PARA *par ) {

    char command[255];
    sprintf ( command, "../msc/initial_data/run_average %d %lf %lf %s %s", par->N, par->xMin, par->xRange, par->project, par->shape );
    if ( !system ( command ) )
        printf ( "# SUCCESS # - created initial data.\n" );

}

static void set_initial_data_burger ( double **Avg, PARA *par ) {

    ASSERT ( Avg );
    ASSERT ( par );

    int i = 0;
    int check;
    double temp;
    double unnecessary;
    char filename[255];

// first cell averages
    sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", par->project, par->shape, par->N );
    FILE *data = fopen ( filename, "r" );   // open file for reading
    ASSERT ( data );

    OUTPUT_EVO ( "opened cell average data" )

    while ( !feof( data ) ) {	// go through file line by line
        check = fscanf ( data, "%lf %lf", &unnecessary, &temp);
        ASSERT ( check );
        Avg[0][i++] = temp;		// write result in array and increment i
    }
    fclose ( data );			// close file

}

void Burger_Init_Data	( FV_DATA *V, PARA *par ) {

    ASSERT ( V );
    ASSERT ( par );

    create_initial_data_burger( par );
    set_initial_data_burger ( V->Avg, par );

    printf ( "# SUCCESS # - initialized wave FV_DATA.\n" );

}

// ### ### ### ### ### ### ###

void Burger_Output_Averages ( double **V, PARA *par, char *name ) {

    int i;
    char filename[255];

    sprintf ( filename, "../data/%s/%s/%s/%d_comp0_%s.txt", par->project, par->shape, par->scheme, par->N, name );

    FILE *fp = fopen ( filename, "a" );

    fprintf ( fp, "\"Time=%lf\n", par->time );

    for ( i = 0; i < par->N; ++i )
        fprintf ( fp, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, V[0][i] );

    fclose ( fp );

}

// ### ### ### ### ### ### ###

void Burger_Output_Averages_single_time ( double **V, PARA *par, char *name ) {

    int i;
    char filename[255];

    sprintf ( filename, "../data/%s/%s/%s/%d_%s_t%.5f.csv", par->project, par->shape, par->scheme, par->N, name, par->time );

    FILE *fp = fopen ( filename, "w" );

    for ( i = 0; i < par->N; ++i )
        fprintf ( fp, "%.15f, %.15f, %.15f\n", par->xMin + (i+0.5)*par->dx, par->time, V[0][i] );

    fclose ( fp );

}

// #### #### ####

void Burger_Flux_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    // first component
    gsl_vector_memcpy ( A->F[0][i], A->Vold[0][i] );
    gsl_vector_mul    ( A->F[0][i], A->Vold[0][i] );
    gsl_vector_scale  ( A->F[0][i], 0.5 );

    par->evalTempF += par->size;
}

void Burger_Flux_Coeff_all ( ADER_DATA *A, PARA *par ) {

    int i;
    for ( i = 0; i < par->N; ++i ) {
        // first component
        gsl_vector_memcpy ( A->F[0][i], A->Vold[0][i] );
        gsl_vector_mul    ( A->F[0][i], A->Vold[0][i] );
        gsl_vector_scale  ( A->F[0][i], 0.5 );

        par->evalTempF += par->size;
    }

}

// #### #### ####

void Burger_Source_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {
// Todo: switch automatically between different source terms in [Montecinos 2015]!?

    /// linear
//    gsl_vector_memcpy   ( A->S[0][i], A->Vold[0][i] );
//    gsl_vector_scale    ( A->S[0][i], -2.0 );

    /// quadratic
//    gsl_vector_memcpy   ( A->S[0][i], A->Vold[0][i] );
//    gsl_vector_mul      ( A->S[0][i], A->S[0][i] );
//    gsl_vector_scale    ( A->S[0][i], -2.0 );

    /// exponential
    int j;
    for ( j = 0; j < par->size; ++j )
        gsl_vector_set ( A->S[0][i], j, exp( -1.0 * gsl_vector_get(A->Vold[0][i], j) ) );

    par->evalTempS += par->size;

}

void Burger_Source_Coeff_all ( ADER_DATA *A, PARA *par ) {
// Todo: switch automatically between different source terms in [Montecinos 2015]!?

    int i;
    int j;
    for ( i = 0; i < par->N; ++i ) {

        /// linear
//      gsl_vector_memcpy   ( A->S[0][i], A->Vold[0][i] );
//      gsl_vector_scale    ( A->S[0][i], -2.0 );

        /// quadratic
//      gsl_vector_memcpy   ( A->S[0][i], A->Vold[0][i] );
//      gsl_vector_mul      ( A->S[0][i], A->S[0][i] );
//      gsl_vector_scale    ( A->S[0][i], -2.0 );

        /// exponential
        for ( j = 0; j < par->size; ++j )
            gsl_vector_set ( A->S[0][i], j, exp( -1.0 * gsl_vector_get(A->Vold[0][i], j) ) );

        par->evalTempS += par->size;

    }

}

// #### #### ####

double Burger_Calc_Speed ( double **V, PARA *par, int i ) {
/// only current cell
    return fabs(V[0][i]);

}

// ### #### ###

void Burger_Convergence_Output ( FV_DATA *V, PARA *par, char *name ) {

    ASSERT ( V );
    ASSERT ( par );
    ASSERT ( name );

    int i;
    char file[255];

    sprintf ( file, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/%s.txt", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name );

    FILE *fp = fopen ( file, "w" );

    for ( i = 0; i < par->N-1; ++i )
        fprintf ( fp, "%.15f %.15f\n", (i+0.5)*par->dx, V->Avg[0][i] );

    fprintf ( fp, "%.15f %.15f", (par->N-0.5)*par->dx, V->Avg[0][(par->N-1)] );

    fclose ( fp );

    printf ( "# SUCCESS # - output convergence file: ../%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/%s.txt\n", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name);

}

// #### #### ####

void Burger_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t ) {

    gsl_vector_memcpy   ( F[0], q[0][i][t] );
    gsl_vector_mul      ( F[0], F[0] );
    gsl_vector_scale    ( F[0], 0.5 );

    par->evalTempF += par->Mp1;

}

void Burger_S_RK ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM ) {
// Todo: switch automatically between different source terms in [Montecinos 2015]!?

    /// linear
//    gsl_vector_memcpy   ( S[0], q[0][i][t] );
//    gsl_vector_scale    ( S[0], -2.0 );

    /// quadratic
//    gsl_vector_memcpy   ( S[0], q[0][i][t] );
//    gsl_vector_mul      ( S[0], S[0] );
//    gsl_vector_scale    ( S[0], -2.0 );

    /// exponential
    int j;
    for ( j = 0; j < par->Mp1; ++j )
        gsl_vector_set ( S[0], j, exp( -1.0 * gsl_vector_get(q[0][i][t], j) ) );

    par->evalTempS += par->Mp1;

}









































