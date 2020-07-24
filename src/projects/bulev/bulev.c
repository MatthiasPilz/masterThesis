#include "ader.h"
#include "bulev.h"

void BuLev_Var_Num ( PARA *par ) {

    par->var = 1;
    printf ( "bulev equation needs %d variables\n", par->var );

}

// ### ### ### ### ### ### ###

static void create_initial_data_bulev ( PARA *par ) {

    char command[255];
    sprintf ( command, "../msc/initial_data/run_average %d %lf %lf %s %s", par->N, par->xMin, par->xRange, par->project, par->shape );
    if ( !system ( command ) )
        printf ( "# SUCCESS # - created initial data.\n" );

}

static void set_initial_data_bulev ( double **Avg, PARA *par ) {

    ASSERT ( Avg );
    ASSERT ( par );

    int i = 0;
    int check;
    double temp;
    char filename[255];

// first cell averages
    sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", par->project, par->shape, par->N );
    FILE *data = fopen ( filename, "r" );   // open file for reading
    ASSERT ( data );

    OUTPUT_EVO ( "opened cell average data" )

    while ( !feof( data ) ) {	// go through file line by line
        check = fscanf ( data, "%lf", &temp);
        ASSERT ( check );
        Avg[0][i++] = temp;		// write result in array and increment i
    }
    fclose ( data );			// close file

}

void BuLev_Init_Data	( FV_DATA *V, PARA *par ) {

    ASSERT ( V );
    ASSERT ( par );

    create_initial_data_bulev( par );
    set_initial_data_bulev ( V->Avg, par );

    printf ( "# SUCCESS # - initialized wave FV_DATA.\n" );

}

// #### #### ####

void BuLev_Flux_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    // need additional helping vector
    gsl_vector *help = gsl_vector_calloc ( par->size );
    gsl_vector *temp = gsl_vector_calloc ( par->size );

    // first component
    gsl_vector_memcpy   ( A->F[0][i], A->Vold[0][i] ); // F = u
    gsl_vector_mul      ( A->F[0][i], A->Vold[0][i] ); // F = u^2
    gsl_vector_scale    ( A->F[0][i], 4.0 ); // F = 4 u^2

    gsl_vector_set_all  ( help, 1.0 ); // help = 1
    gsl_vector_sub      ( help, A->Vold[0][i] ); // help = 1 - u

    gsl_vector_memcpy   ( temp, help ); // temp = 1 - u

    gsl_vector_mul      ( help, temp ); // FEHLERQUELLE!? help = (1-u)^2
    gsl_vector_add      ( help, A->F[0][i] ); // help = 4u^2 + (1-u)^2

    gsl_vector_div      ( A->F[0][i], help );

}

void BuLev_Flux_Coeff_all ( ADER_DATA *A, PARA *par ) {

    int i;
    // need additional helping vector
    gsl_vector *help = gsl_vector_calloc ( par->size );
    gsl_vector *temp = gsl_vector_calloc ( par->size );

    for ( i = 0; i < par->N; ++i ) {

    // first component
        gsl_vector_memcpy   ( A->F[0][i], A->Vold[0][i] ); // F = u
        gsl_vector_mul      ( A->F[0][i], A->Vold[0][i] ); // F = u^2
        gsl_vector_scale    ( A->F[0][i], 4.0 ); // F = 4 u^2

        gsl_vector_set_all  ( help, 1.0 ); // help = 1
        gsl_vector_sub      ( help, A->Vold[0][i] ); // help = 1 - u

        gsl_vector_memcpy   ( temp, help ); // temp = 1 - u

        gsl_vector_mul      ( help, temp ); // FEHLERQUELLE!? help = (1-u)^2
        gsl_vector_add      ( help, A->F[0][i] ); // help = 4u^2 + (1-u)^2

        gsl_vector_div      ( A->F[0][i], help );
    }

}

// #### #### ####

double BuLev_Calc_Speed ( double **V, PARA *par, int i ) {

    double result;

    result = (8 * V[0][i] * ( 1-V[0][i] )) / ( pow( 1 - 2*V[0][i] + 5*V[0][i]*V[0][i],2 ) );

    return result;

}









































