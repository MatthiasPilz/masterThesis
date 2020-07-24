#include "ader.h"
#include "euler.h"

void Euler_Var_Num ( PARA *par ) {

    par->var = 3;
    printf ( "euler equations needs %d variables\n", par->var );

}

// ### ### ### ### ### ### ###

static void create_initial_data_euler ( PARA *par ) {

    char command[255];
    sprintf ( command, "../msc/initial_data/run_average %d %lf %lf %s %s", par->N, par->xMin, par->xRange, par->project, par->shape );
    printf ( "%s\n", command );

    if ( !system ( command ) )
        printf ( "# SUCCESS # - created initial data.\n" );

}

static void set_initial_data_euler ( double **Avg, PARA *par ) {

    ASSERT ( Avg );
    ASSERT ( par );

    int i = 0;
    int check;
    double temp;
    double unnecessary;
    char filename[255];

    // helping arrays...
    double *R = calloc ( par->N, sizeof ( double ) );
    double *V = calloc ( par->N, sizeof ( double ) );
    double *P = calloc ( par->N, sizeof ( double ) );

// first read primitive cell averages
    /// density
    sprintf ( filename, "../math/initial_data/%s/%s/%d_density.txt", par->project, par->shape, par->N );
    FILE *density = fopen ( filename, "r" );   // open file for reading
    ASSERT ( density );

    OUTPUT_EVO ( "opened cell average data" )

    while ( !feof( density ) ) {	// go through file line by line
        check = fscanf ( density, "%lf %lf", &unnecessary, &temp);
        ASSERT ( check );
        R[i++] = temp;		// write result in array and increment i
    }
    fclose ( density );			// close file

    /// pressure
    i = 0;
    sprintf ( filename, "../math/initial_data/%s/%s/%d_pressure.txt", par->project, par->shape, par->N );
    FILE *pressure = fopen ( filename, "r" );   // open file for reading
    ASSERT ( pressure );

    OUTPUT_EVO ( "opened cell average data" )

    while ( !feof( pressure ) ) {	// go through file line by line
        check = fscanf ( pressure, "%lf %lf", &unnecessary, &temp);
        ASSERT ( check );
        P[i++] = temp;		// write result in array and increment i
    }
    fclose ( pressure );			// close file

    /// velocity
    i = 0;
    sprintf ( filename, "../math/initial_data/%s/%s/%d_velocity.txt", par->project, par->shape, par->N );
    FILE *velocity = fopen ( filename, "r" );   // open file for reading
    ASSERT ( velocity );

    OUTPUT_EVO ( "opened cell average data" )

    while ( !feof( velocity ) ) {	// go through file line by line
        check = fscanf ( velocity, "%lf %lf", &unnecessary, &temp);
        ASSERT ( check );
        V[i++] = temp;		// write result in array and increment i
    }
    fclose ( velocity );			// close file


// set conserved variables
    for ( i = 0; i < par->N; ++i ) {
        Avg[0][i] = R[i];
        Avg[1][i] = R[i]*V[i];
        Avg[2][i] = R[i] * ( 0.5*V[i]*V[i] + P[i]/((par->e_gamma-1.0) * R[i]) );
    }

    free ( R );
    free ( P );
    free ( V );

}

void Euler_Init_Data ( FV_DATA *V, PARA *par ) {

    ASSERT ( V );
    ASSERT ( par );

    create_initial_data_euler       ( par );
    set_initial_data_euler          ( V->Avg, par );

    printf ( "# SUCCESS # - initialized wave FV_DATA.\n" );

}

// ### ### ### ### ### ### ###

void Euler_Output_Averages ( double **V, PARA *par, char *name ) {

    int i;
    char filename[255];

// 1st component
    sprintf ( filename, "../data/%s/%s/%s/%d_comp0.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp = fopen ( filename, "a" );
    fprintf ( fp, "\"Time=%lf\n", par->time );
    for ( i = 0; i < par->N; ++i )
        fprintf ( fp, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, V[0][i] );

    fclose ( fp );

// 2nd component
    sprintf ( filename, "../data/%s/%s/%s/%d_comp1.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp1 = fopen ( filename, "a" );
    fprintf ( fp1, "\"Time=%lf\n", par->time );
    for ( i = 0; i < par->N; ++i )
        fprintf ( fp1, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, V[1][i]/V[0][i] );

    fclose ( fp1 );

// 3rd component
    sprintf ( filename, "../data/%s/%s/%s/%d_comp2.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp2 = fopen ( filename, "a" );
    fprintf ( fp2, "\"Time=%lf\n", par->time );
    for ( i = 0; i < par->N; ++i )
        fprintf ( fp2, "%.15f %.15f\n", par->xMin + (i+0.5)*par->dx, ( par->e_gamma - 1.0 ) * ( V[2][i] - 0.5 * ( V[1][i] * V[1][i] / V[0][i] ) ) );

    fclose ( fp2 );

}

// ### ### ### ### ### ### ###

void Euler_Flux_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    gsl_vector *help = gsl_vector_calloc ( par->size );

    // first component
    gsl_vector_memcpy ( A->F[0][i], A->Vold[1][i] );

    // second component
    gsl_vector_memcpy   ( help, A->Vold[1][i] );
    gsl_vector_mul      ( help, A->Vold[1][i] );
    gsl_vector_div      ( help, A->Vold[0][i] );
    gsl_vector_scale    ( help, 0.5*( 3.0 - par->e_gamma ) );

    gsl_vector_memcpy   ( A->F[1][i], A->Vold[2][i] );
    gsl_vector_scale    ( A->F[1][i], (par->e_gamma - 1.0) );

    gsl_vector_add      ( A->F[1][i], help );

    // third component
    gsl_vector_set_zero ( help );
    gsl_vector_memcpy   ( help, A->Vold[1][i] );
    gsl_vector_mul      ( help, A->Vold[2][i] );
    gsl_vector_div      ( help, A->Vold[0][i] );
    gsl_vector_scale    ( help, par->e_gamma );

    gsl_vector_memcpy   ( A->F[2][i], A->Vold[1][i] );
    gsl_vector_mul      ( A->F[2][i], A->Vold[1][i] );
    gsl_vector_mul      ( A->F[2][i], A->Vold[1][i] );
    gsl_vector_div      ( A->F[2][i], A->Vold[0][i] );
    gsl_vector_div      ( A->F[2][i], A->Vold[0][i] );
    gsl_vector_scale    ( A->F[2][i], -0.5*( par->e_gamma - 1.0 ) );

    gsl_vector_add      ( A->F[2][i], help );

    par->evalTempF += par->var * par->size;

    gsl_vector_free ( help );
}

void Euler_Flux_Coeff_all ( ADER_DATA *A, PARA *par ) {

    int i;
    gsl_vector *help = gsl_vector_calloc ( par->size );

    for ( i = 0; i < par->N; ++i ) {
        gsl_vector_set_zero ( help );

    // first component
        gsl_vector_memcpy ( A->F[0][i], A->Vold[1][i] );

    // second component
        gsl_vector_memcpy   ( help, A->Vold[1][i] );
        gsl_vector_mul      ( help, A->Vold[1][i] );
        gsl_vector_div      ( help, A->Vold[0][i] );
        gsl_vector_scale    ( help, 0.5*( 3.0 - par->e_gamma ) );

        gsl_vector_memcpy   ( A->F[1][i], A->Vold[2][i] );
        gsl_vector_scale    ( A->F[1][i], (par->e_gamma - 1.0) );

        gsl_vector_add      ( A->F[1][i], help );

    // third component
        gsl_vector_set_zero ( help );
        gsl_vector_memcpy   ( help, A->Vold[1][i] );
        gsl_vector_mul      ( help, A->Vold[2][i] );
        gsl_vector_div      ( help, A->Vold[0][i] );
        gsl_vector_scale    ( help, par->e_gamma );

        gsl_vector_memcpy   ( A->F[2][i], A->Vold[1][i] );
        gsl_vector_mul      ( A->F[2][i], A->Vold[1][i] );
        gsl_vector_mul      ( A->F[2][i], A->Vold[1][i] );
        gsl_vector_div      ( A->F[2][i], A->Vold[0][i] );
        gsl_vector_div      ( A->F[2][i], A->Vold[0][i] );
        gsl_vector_scale    ( A->F[2][i], -0.5*( par->e_gamma - 1.0 ) );

        gsl_vector_add      ( A->F[2][i], help );

        par->evalTempF += par->var * par->size;

    }

    gsl_vector_free ( help );
}

// ### ### ### ### ### ### ###

double Euler_Calc_Speed_av ( double **V, PARA *par, int i ) {

    ASSERT ( V );

    double s1;
    double s2;
    double s3;

    double temp1;
    double temp2;
    double temp3;

    if ( i != 0 && i != par->N-1 ) {
        temp1 = ( par->e_gamma - 1.0 ) * ( V[2][i] - 0.5 * ( V[1][i] * V[1][i] / V[0][i] ) );
        temp2 = ( par->e_gamma - 1.0 ) * ( V[2][i+1] - 0.5 * ( V[1][i+1] * V[1][i+1] / V[0][i+1] ) );
        temp3 = ( par->e_gamma - 1.0 ) * ( V[2][i-1] - 0.5 * ( V[1][i-1] * V[1][i-1] / V[0][i-1] ) );

        s1 = ( (V[1][ i ]/V[0][ i ]) + sqrt ( par->e_gamma * temp1 / V[0][ i ] ) );
        s2 = ( (V[1][i+1]/V[0][i+1]) + sqrt ( par->e_gamma * temp2 / V[0][i+1] ) );
        s3 = ( (V[1][i-1]/V[0][i-1]) + sqrt ( par->e_gamma * temp3 / V[0][i-1] ) );

    }
    else if ( strcmp( par->boundary, "periodic" ) == 0 ) {
        if ( i == 0 ) {
            temp1 = ( par->e_gamma - 1.0 ) * ( V[2][ i ] - 0.5 * ( V[1][ i ] * V[1][ i ] / V[0][ i ] ) );
            temp2 = ( par->e_gamma - 1.0 ) * ( V[2][i+1] - 0.5 * ( V[1][i+1] * V[1][i+1] / V[0][i+1] ) );
            temp3 = ( par->e_gamma - 1.0 ) * ( V[2][par->N-1] - 0.5 * ( V[1][par->N-1] * V[1][par->N-1] / V[0][par->N-1] ) );

            s1 = ( (V[1][ i ]/V[0][ i ]) + sqrt ( par->e_gamma * temp1 / V[0][ i ] ) );
            s2 = ( (V[1][i+1]/V[0][i+1]) + sqrt ( par->e_gamma * temp2 / V[0][i+1] ) );
            s3 = ( (V[1][par->N-1]/V[0][par->N-1]) + sqrt ( par->e_gamma * temp3 / V[0][par->N-1] ) );
        }
        else if ( i == par->N-1 ) {
            temp1 = ( par->e_gamma - 1.0 ) * ( V[2][ i ] - 0.5 * ( V[1][ i ] * V[1][ i ] / V[0][ i ] ) );
            temp2 = ( par->e_gamma - 1.0 ) * ( V[2][ 0 ] - 0.5 * ( V[1][ 0 ] * V[1][ 0 ] / V[0][ 0 ] ) );
            temp3 = ( par->e_gamma - 1.0 ) * ( V[2][i-1] - 0.5 * ( V[1][i-1] * V[1][i-1] / V[0][i-1] ) );

            s1 = ( (V[1][ i ]/V[0][ i ]) + sqrt ( par->e_gamma * temp1 / V[0][ i ] ) );
            s2 = ( (V[1][ 0 ]/V[0][ 0 ]) + sqrt ( par->e_gamma * temp2 / V[0][ 0 ] ) );
            s3 = ( (V[1][i-1]/V[0][i-1]) + sqrt ( par->e_gamma * temp3 / V[0][i-1] ) );
        }
    }
    else if ( strcmp( par->boundary, "none" ) == 0 ) {
        if ( i == 0 ) {
            temp1 = ( par->e_gamma - 1.0 ) * ( V[2][0] - 0.5 * ( V[1][0] * V[1][0] / V[0][0] ) );
            s1 = ( (V[1][ 0 ]/V[0][ 0 ]) + sqrt ( par->e_gamma * temp1 / V[0][ 0 ] ) );
            s2 = s1;
            s3 = s1;
        }
        else if ( i == par->N-1 ) {
            temp1 = ( par->e_gamma - 1.0 ) * ( V[2][par->N-1] - 0.5 * ( V[1][par->N-1] * V[1][par->N-1] / V[0][par->N-1] ) );
            s1 = ( (V[1][ par->N-1 ]/V[0][ par->N-1 ]) + sqrt ( par->e_gamma * temp1 / V[0][ par->N-1 ] ) );
            s2 = s1;
            s3 = s1;
        }
    }
    else
        printf ( "*** ERROR - Euler speed calculation\n" );

    return fabs( (s1 + s2 + s3) / 3.0 );
}

//### ### ###

void Euler_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t ) {

        gsl_vector *help = gsl_vector_calloc ( par->Mp1 );
        gsl_vector_set_zero ( help );

    // first component
        gsl_vector_memcpy ( F[0], q[1][i][t] );

    // second component
        gsl_vector_memcpy   ( help, q[1][i][t] );
        gsl_vector_mul      ( help, q[1][i][t] );
        gsl_vector_div      ( help, q[0][i][t] );
        gsl_vector_scale    ( help, 0.5*( 3.0 - par->e_gamma ) );

        gsl_vector_memcpy   ( F[1], q[2][i][t] );
        gsl_vector_scale    ( F[1], (par->e_gamma - 1.0) );

        gsl_vector_add      ( F[1], help );

    // third component
        gsl_vector_set_zero ( help );
        gsl_vector_memcpy   ( help, q[1][i][t] );
        gsl_vector_mul      ( help, q[2][i][t] );
        gsl_vector_div      ( help, q[0][i][t] );
        gsl_vector_scale    ( help, par->e_gamma );

        gsl_vector_memcpy   ( F[2], q[1][i][t] );
        gsl_vector_mul      ( F[2], q[1][i][t] );
        gsl_vector_mul      ( F[2], q[1][i][t] );
        gsl_vector_div      ( F[2], q[0][i][t] );
        gsl_vector_div      ( F[2], q[0][i][t] );
        gsl_vector_scale    ( F[2], -0.5*( par->e_gamma - 1.0 ) );

        gsl_vector_add      ( F[2], help );

        par->evalTempF += par->var * par->Mp1;

        gsl_vector_free ( help );
}

// ### #### ###

void Euler_Convergence_Output ( FV_DATA *V, PARA *par, char *name ) {

    ASSERT ( V );
    ASSERT ( par );
    ASSERT ( name );

    int i;
    char file[255];

    sprintf ( file, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/%s.txt", par->project, par->M, par->shape, par->scheme, par->CFL, par->N, name );

    FILE *fp = fopen ( file, "w" );

    for ( i = 0; i < par->N-1; ++i )
        fprintf ( fp, "%.15f %.15f\n", (i+0.5)*par->dx, V->Avg[0][i] );

    fprintf ( fp, "%.15f %.15f", (par->N-0.5)*par->dx, V->Avg[0][(par->N-1)] );

    fclose ( fp );

    printf ( "# SUCCESS # - output convergence file: ../%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/%s.txt\n", par->project, par->M, par->shape, par->scheme, par->CFL, par->N, name);

}







