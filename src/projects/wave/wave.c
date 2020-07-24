#include "ader.h"
#include "wave.h"

void Wave_Var_Num ( PARA *par ) {

	par->var = 3;
	printf ( "wave equation needs %d variables\n", par->var );

}

// ### ### ### ### ### ### ###

static void create_initial_data_wave ( PARA *par ) {

    char command[255];
    sprintf ( command, "../msc/initial_data/run_average %d %lf %lf %s %s", par->N, par->xMin, par->xRange, par->project, par->shape );
    if ( !system ( command ) )
        printf ( "# SUCCESS # - created initial data.\n" );

}

static void set_initial_data_wave ( double **Avg, PARA *par ) {

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


// now spatial derivative
    i = 0;
    sprintf ( filename, "../math/initial_data/%s/%s/%d_ICD.txt", par->project, par->shape, par->N );
    FILE *fp = fopen ( filename, "r" );   // open file for reading
    ASSERT ( data );

    OUTPUT_EVO ( "opened derivative data" )

    while ( !feof( fp ) ) {	// go through file line by line
        check = fscanf ( fp, "%lf %lf", &unnecessary, &temp);
        ASSERT ( check );
        Avg[1][i++] = temp;			// write result in array and increment i
    }
    fclose ( fp );			// close file

}

void Wave_Init_Data	( FV_DATA *V, PARA *par ) {

    ASSERT ( V );
    ASSERT ( par );

    create_initial_data_wave( par );
    set_initial_data_wave ( V->Avg, par );

    printf ( "# SUCCESS # - initialized wave FV_DATA.\n" );

}

// ### ### ### ### ### ### ###

void Wave_Output_Averages ( double **V, PARA *par, char *name ) {

    int i, j;
    char filename[255];

    for ( j = 0; j < par->var; ++j ) {
        sprintf ( filename, "../data/%s/%s/%s/%d_comp%d_%s.txt", par->project, par->shape, par->scheme, par->N, j, name );

        FILE *fp = fopen ( filename, "a" );

        fprintf ( fp, "\"Time=%lf\n", par->time );

        for ( i = 0; i < par->N; ++i )
            fprintf ( fp, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, V[j][i] );

        fclose ( fp );
    }

}

// #### #### ####

void Wave_Flux_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    // first component
    gsl_vector_set_zero ( A->F[0][i] );

    // second component
    gsl_vector_memcpy ( A->F[1][i], A->Vold[2][i] );
    gsl_vector_scale  ( A->F[1][i], -1.0 );

    // third component
    gsl_vector_memcpy ( A->F[2][i], A->Vold[1][i] );
    gsl_vector_scale  ( A->F[2][i], -par->speed*par->speed );

    par->evalTempF += par->var*par->size;

}

void Wave_Flux_Coeff_all ( ADER_DATA *A, PARA *par ) {

    int i;
    for ( i = 0; i < par->N; ++i ) {

    // first component
        gsl_vector_set_zero ( A->F[0][i] );

    // second component
        gsl_vector_memcpy ( A->F[1][i], A->Vold[2][i] );
        gsl_vector_scale  ( A->F[1][i], -1.0 );

    // third component
        gsl_vector_memcpy ( A->F[2][i], A->Vold[1][i] );
        gsl_vector_scale  ( A->F[2][i], -par->speed*par->speed );

        par->evalTempF += par->var*par->size;

    }

}

// #### #### ####

void Wave_Source_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    // first component
    gsl_vector_memcpy ( A->S[0][i], A->Vold[2][i] );

    // second and third
    gsl_vector_set_zero ( A->S[1][i] );
    gsl_vector_set_zero ( A->S[2][i] );

    par->evalTempS += par->var*par->size;

}

void Wave_Source_Coeff_all ( ADER_DATA *A, PARA *par ) {

    int i;
    for ( i = 0; i < par->N; ++i ) {

    // first component
        gsl_vector_memcpy ( A->S[0][i], A->Vold[2][i] );

    // second and third
        gsl_vector_set_zero ( A->S[1][i] );
        gsl_vector_set_zero ( A->S[2][i] );

        par->evalTempS += par->var*par->size;

    }

}

// #### #### ####

double Wave_Calc_Speed ( double **V, PARA *par, int i ) {

    return par->speed;

}

// ### #### ###

void Wave_Convergence_Output ( FV_DATA *V, PARA *par, char *name ) {

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

// ### ### ### ### ### ### ### ###

void Wave_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t ) {

    // first component
    gsl_vector_set_zero ( F[0] );

    // second component
    gsl_vector_memcpy ( F[1], q[2][i][t] );
    gsl_vector_scale  ( F[1], -1.0 );

    // third component
    gsl_vector_memcpy ( F[2], q[1][i][t] );
    gsl_vector_scale  ( F[2], -par->speed*par->speed );

    par->evalTempF += par->var*par->Mp1;

}

void Wave_S_RK ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM ) {


    // first component
    gsl_vector_memcpy ( S[0], q[2][i][t] );

    // second and third
    gsl_vector_set_zero ( S[1] );
    gsl_vector_set_zero ( S[2] );

    par->evalTempS += par->var*par->Mp1;

}









































