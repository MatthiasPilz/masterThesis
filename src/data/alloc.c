#include "ader.h"
#include "data.h"


void Allocate_FV_DATA_general ( FV_DATA *V, PARA *par ) {
/*
 * V    - FV_DATA
 * par  - parameter
 *
 * allocate memory for all arrays and matrices in FV_DATA
 */

    ASSERT ( par );
    ASSERT ( V );

    int k;

    // cell average data arrays
    V->Avg = malloc ( par->var * sizeof ( &V->Avg ) );
    for ( k = 0; k < par->var; ++k )
        V->Avg[k] = calloc ( par->N, sizeof ( double ) );

    V->Avg_adm = malloc ( par->var * sizeof ( &V->Avg_adm ) );
    for ( k = 0; k < par->var; ++k )
        V->Avg_adm[k] = calloc ( par->N, sizeof ( double ) );

    V->Flux = malloc ( par->var * sizeof ( &V->Flux ) );
    for ( k = 0; k < par->var; ++k )
        V->Flux[k] = calloc ( par->N+1, sizeof ( double ) );

    if ( par->source_flag ) {
        V->Source = malloc ( par->var * sizeof ( &V->Source ) );
        for ( k = 0; k < par->var; ++k )
            V->Source[k] = calloc ( par->N, sizeof ( double ) );
    }

    if ( strcmp( par->boundary, "hack") == 0 ) {
        V->hack_Flux = malloc ( par->var * sizeof ( &V->hack_Flux ) );
        for ( k = 0; k < par->var; ++k )
            V->hack_Flux[k] = calloc ( 1, sizeof ( double ) );
    }

    // precomputed matrices
    // time averaged flux matrix eq.(4) in [1]
    V->F1_1 = gsl_vector_calloc ( par->size );
    V->F1_0 = gsl_vector_calloc ( par->size );

    // time and space averaged source matrix eq.(7) in [1]
    if ( par->source_flag )
        V->S1   = gsl_vector_calloc ( par->size );

    printf ( "# SUCCESS # - allocated memory for FV_DATA struct.\n" );

}

// ### ### ### ### ### ### ###

void Allocate_WENO_DATA_general ( WENO_DATA *W, PARA *par ) {
/*
 * W    - WENO_DATA
 * par  - parameter
 *
 * allocate all arrays, vectors and matrices in WENO_DATA
 */

    ASSERT ( W );
    ASSERT ( par );
    ASSERT ( par->var > 0 );

    int i, j, k;

//  gslData
    W->gslData = malloc ( par->var * sizeof ( &W->gslData ) );
    for ( k = 0; k < par->var; ++k )
        W->gslData[k] = malloc ( par->N * sizeof ( &W->gslData[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            W->gslData[k][i] = malloc ( par->Mp1 * sizeof ( &W->gslData[k][i] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            for ( j = 0; j < par->Mp1; ++j )
                W->gslData[k][i][j] = gsl_vector_calloc ( par->Mp1 );

// gslWs
    W->gslWs = malloc ( par->var * sizeof ( &W->gslWs ) );
    for ( k = 0; k < par->var; ++k )
        W->gslWs[k] = malloc ( par->N * sizeof ( &W->gslWs[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            W->gslWs[k][i] = malloc ( par->Mp1 * sizeof ( &W->gslWs[k][i] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            for ( j = 0; j < par->Mp1; ++j )
                W->gslWs[k][i][j] = gsl_vector_calloc ( par->Mp1 );

// reconstruction matrices
    W->Rec = malloc ( par->Mp1 * sizeof ( &W->Rec ) );
    for ( j = 0; j < par->n_stencil; ++j )
        W->Rec[j] = gsl_matrix_calloc ( par->Mp1, par->Mp1 );

// l_s
    W->l_s = malloc ( par->n_stencil * sizeof ( double ) );

// shift
    W->shift = calloc ( par->n_stencil, sizeof ( int ) );

// OIM
    W->OIM = malloc ( par->Mp1 * sizeof ( &W->OIM ) );
    for ( j = 0; j < par->Mp1; ++j )
        W->OIM[j] = calloc ( par->Mp1, sizeof ( double ) );

// w_s
    W->w_s = malloc ( par->N * sizeof ( &W->w_s ) );
    for ( i = 0; i < par->N; ++i )
        W->w_s[i] = calloc ( par->Mp1, sizeof ( double ) );

// w_tilde
    W->w_tilde = malloc ( par->N * sizeof ( &W->w_tilde ) );
    for ( i = 0; i < par->N; ++i )
        W->w_tilde[i] = calloc ( par->Mp1, sizeof ( double ) );

// sigma
    W->sigma = malloc ( par->N * sizeof ( &W->sigma ) );
    for ( i = 0; i < par->N; ++i )
        W->sigma[i] = calloc ( par->Mp1, sizeof ( double ) );

// WENO coefficients
    W->WT = malloc ( par->var * sizeof ( &W->WT ) );
    for ( k = 0; k < par->var; ++k )
        W->WT[k] = malloc ( par->N * sizeof ( &W->WT[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            W->WT[k][i] = gsl_vector_calloc( par->Mp1 );

    printf ( "# SUCCESS # - allocated memory for WENO_DATA struct.\n" );

}

// ### ### ### ### ### ### ###

void Allocate_ADER_DATA_general ( ADER_DATA *A, PARA *par ) {
/*
 * A    - ADER_DATA
 * par  - parameter
 *
 * allocate all vectors and matrices in ADER_DATA
 */

    ASSERT ( A );
    ASSERT ( par );
    ASSERT ( par->var > 0 );

    int i, k, var;

// Vold
    A->Vold = malloc ( par->var * sizeof ( &A->Vold ) );
    for ( k = 0; k < par->var; ++k )
        A->Vold[k] = malloc ( par->N * sizeof ( &A->Vold[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            A->Vold[k][i] = gsl_vector_calloc ( par->size );

// Vnew
    A->Vnew = malloc ( par->var * sizeof ( &A->Vnew ) );
    for ( k = 0; k < par->var; ++k )
        A->Vnew[k] = malloc ( par->N * sizeof ( &A->Vnew[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            A->Vnew[k][i] = gsl_vector_calloc ( par->size );

// F - flux
    A->F = malloc ( par->var * sizeof ( &A->F ) );
    for ( k = 0; k < par->var; ++k )
        A->F[k] = malloc ( par->N * sizeof ( &A->F[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            A->F[k][i] = gsl_vector_calloc ( par->size );

// hack_F - flux for inside of TOV star
    if ( strcmp( par->boundary, "hack") == 0 ) {
        A->hack_F = malloc ( par->var * sizeof ( &A->hack_F ) );
        for ( k = 0; k < par->var; ++k )
            A->hack_F[k] = malloc ( 1 * sizeof ( &A->hack_F[k] ) );

        for ( k = 0; k < par->var; ++k )
            A->hack_F[k][0] = gsl_vector_calloc ( par->size );
    }

// S - source
    A->S = malloc ( par->var * sizeof ( &A->S ) );
    for ( k = 0; k < par->var; ++k )
        A->S[k] = malloc ( par->N * sizeof ( &A->S[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            A->S[k][i] = gsl_vector_calloc ( par->size );

// matrix I
    A->I = gsl_matrix_calloc ( par->size, par->Mp1 );

// matrix B
    A->B = gsl_matrix_calloc ( par->size, par->size );

// matrix MS
    A->Ms = gsl_matrix_calloc( par->size, par->size );

// LAM
    A->LAM = calloc ( par->Mp1, sizeof ( double ) );

// matrix D
    A->D = gsl_matrix_calloc ( par->Mp1, par->Mp1 );


/// RK_stuff
  if ( strcmp ( par->scheme, "RK" ) == 0 ) {

  // Butcher Aij Matrix
        A->Aij = malloc ( par->stages * sizeof ( &A->Aij ) );
        for ( k = 0; k < par->stages; ++k )
            A->Aij[k] = calloc ( par->stages, sizeof ( double ) );


  // k_RK
        A->k_RK = malloc ( par->var * sizeof ( &A->k_RK ) );
        for ( var = 0; var < par->var; ++var )
            A->k_RK[var] = malloc ( par->N * sizeof ( &A->k_RK[var] ) );

        for ( var = 0; var < par->var; ++var )
            for ( i = 0; i < par->N; ++i )
                A->k_RK[var][i] = malloc ( par->stages * sizeof ( &A->k_RK[var][i] ) );

        for ( var = 0; var < par->var; ++var )
            for ( i = 0; i < par->N; ++i )
                for ( k = 0; k < par->stages; ++k )
                    A->k_RK[var][i][k] = gsl_vector_calloc ( par->Mp1 );

  // q_RK
        A->q_RK = malloc ( par->var * sizeof ( &A->q_RK ) );
        for ( var = 0; var < par->var; ++var )
            A->q_RK[var] = malloc ( par->N * sizeof ( &A->q_RK[var] ) );

        for ( var = 0; var < par->var; ++var )
            for ( i = 0; i < par->N; ++i )
                A->q_RK[var][i] = malloc ( par->stages * sizeof ( &A->q_RK[var][i] ) );

        for ( var = 0; var < par->var; ++var )
            for ( i = 0; i < par->N; ++i )
                for ( k = 0; k < par->stages; ++k )
                    A->q_RK[var][i][k] = gsl_vector_calloc ( par->Mp1 );

  // F_RK
        A->F_RK = malloc ( par->var * sizeof ( &A->F_RK ) );
        for ( var = 0; var < par->var; ++var )
            A->F_RK[var] = gsl_vector_calloc ( par->Mp1 );

  // S_RK
        A->S_RK = malloc ( par->var * sizeof ( &A->S_RK ) );
        for ( var = 0; var < par->var; ++var )
            A->S_RK[var] = gsl_vector_calloc ( par->Mp1 );

  // Dense Output
        A->Dense = malloc ( par->Mp1 * sizeof ( &A->Dense ) );
        for ( k = 0; k < par->Mp1; ++k )
            A->Dense[k] = calloc ( par->stages, sizeof ( double ) );

    }



    printf ( "# SUCCESS # - allocated memory for ADER_DATA struct.\n" );

}

// ### ### ### ### ### ### ###

void Allocate_DG_DATA_general ( DG_DATA *DG, PARA *par ) {
/*
 * DG   - DG_DATA
 * par  - parameter
 *
 * allocate all vectors and matrices in DG_DATA
 */

    ASSERT ( DG );
    ASSERT ( par );

    int i, k;

// U
    DG->U = malloc ( par->var * sizeof ( &DG->U ) );
    for ( k = 0; k < par->var; ++k )
       DG->U[k] = malloc ( par->N * sizeof ( &DG->U[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            DG->U[k][i] = gsl_vector_calloc( par->Mp1 );

// pU
    DG->pU = malloc ( par->var * sizeof ( &DG->pU ) );
    for ( k = 0; k < par->var; ++k )
       DG->pU[k] = malloc ( par->N * sizeof ( &DG->pU[k] ) );

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            DG->pU[k][i] = gsl_vector_calloc( par->Mp1 );

// V
    DG->V = malloc ( par->var * sizeof ( &DG->V ) );
    for ( k = 0; k < par->var; ++k )
        DG->V[k] = calloc ( par->N*par->NSC, sizeof ( double ) );

// V_grid
    DG->V_grid = malloc ( par->var * sizeof ( &DG->V_grid ) );
    for ( k = 0; k < par->var; ++k )
        DG->V_grid[k] = calloc ( par->N, sizeof ( double ) );

// V_grid_old
    DG->V_grid_pre = malloc ( par->var * sizeof ( &DG->V_grid_pre ) );
    for ( k = 0; k < par->var; ++k )
        DG->V_grid_pre[k] = calloc ( par->N, sizeof ( double ) );

// beta
    DG->beta = malloc ( par->var * sizeof ( &DG->beta ) );
    for ( k = 0; k < par->var; ++k )
        DG->beta[k] = calloc ( par->N, sizeof ( int ) );

// beta_sub
    DG->beta_sub = malloc ( par->var * sizeof ( &DG->beta_sub ) );
    for ( k = 0; k < par->var; ++k )
        DG->beta_sub[k] = calloc ( par->NSC*par->N, sizeof ( int ) );

// beta_prev
    DG->beta_prev = malloc ( par->var * sizeof ( &DG->beta_prev ) );
    for ( k = 0; k < par->var; ++k )
        DG->beta_prev[k] = calloc ( par->N, sizeof ( int ) );

// Int
    DG->Int = gsl_matrix_calloc ( par->NSC, par->Mp1 );    

// Rec
    DG->Rec = gsl_matrix_calloc ( par->Mp1, par->NSC );

// AVG
    DG->AVG = calloc ( par->Mp1, sizeof ( double ) );

// LAM
    DG->LAM = calloc ( par->Mp1, sizeof ( double ) );

// AVG
    DG->BLow = calloc ( par->Mp1, sizeof ( double ) );

// AVG
    DG->BUp = calloc ( par->Mp1, sizeof ( double ) );

// matrix F
    DG->F   = gsl_matrix_calloc ( par->Mp1, par->size );

// matrix G00
    DG->G00 = gsl_matrix_calloc ( par->Mp1, par->size );

// matrix G10
    DG->G10 = gsl_matrix_calloc ( par->Mp1, par->size );

// matrix G01
    DG->G01 = gsl_matrix_calloc ( par->Mp1, par->size );

// matrix G11
    DG->G11 = gsl_matrix_calloc ( par->Mp1, par->size );

// matrix H
    if ( par->source_flag )
        DG->H   = gsl_matrix_calloc ( par->Mp1, par->size );

    printf ( "# SUCCESS # - allocated memory for DG_DATA struct.\n" );

}




























