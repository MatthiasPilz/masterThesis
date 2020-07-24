#include "ader.h"
#include "data.h"

void Free_FV_DATA ( FV_DATA *V, PARA *par ) {
/*
 * V   - FiniteVolume Data
 * par - parameter
 *
 * free all allocated arrays and gsl_vectors
 */

    ASSERT ( V );
    ASSERT ( par );

    int i;

    // Avg, Flux, Source
    for ( i = 0; i < par->var; ++i )
        free ( V->Avg[i] );

    for ( i = 0; i < par->var; ++i )
        free ( V->Flux[i] );

    if ( par->source_flag ) {
        gsl_vector_free ( V->S1 );

        for ( i = 0; i < par->var; ++i )
            free ( V->Source[i] );
    }

    // gsl_vectors
    gsl_vector_free ( V->F1_0 );
    gsl_vector_free ( V->F1_1 );

    printf ( "# SUCCESS # - freed general FV_DATA.\n" );

}

// ### ### ###

void Free_WENO_DATA ( WENO_DATA *W, PARA *par ) {
/*
 * W   - WENO-Data
 * par - parameter
 *
 * free all allocated arrays and gsl_vectors
 */

    ASSERT ( W );
    ASSERT ( par );
    ASSERT ( par->var > 0 );

    int i, j, k;

    for ( i = 0; i < par->var; ++i )
        for ( j = 0; j < par->N; ++j )
            for ( k = 0; k < par->Mp1; ++k )
                gsl_vector_free ( W->gslData[i][j][k] );

    for ( i = 0; i < par->var; ++i )
        for ( j = 0; j < par->N; ++j )
            for ( k = 0; k < par->Mp1; ++k )
                gsl_vector_free ( W->gslWs[i][j][k] );

    for ( i = 0; i < par->n_stencil; ++i )
        gsl_matrix_free ( W->Rec[i] );

    free ( W->l_s );
    free ( W->shift );

    for ( i = 0; i < par->N; ++i )
        free ( W->w_s[i] );

    for ( i = 0; i < par->N; ++i )
        free ( W->w_tilde[i] );

    for ( i = 0; i < par->N; ++i )
        free ( W->sigma[i] );

    for ( i = 0; i < par->var; ++i )
        for ( j = 0; j < par->N; ++j )
            gsl_vector_free ( W->WT[i][j] );

    printf ( "# SUCCESS # - freed general WENO_DATA.\n" );

}

// ### ### ###

void Free_ADER_DATA ( ADER_DATA *A, PARA *par ) {
/*
 * A   - ADER-Data
 * par - parameter
 *
 * free all allocated arrays and gsl_vectors
 */

    ASSERT ( A );
    ASSERT ( par );

    int i, j;

    for ( i = 0; i < par->var; ++i )
        for ( j = 0; j < par->N; ++j )
            gsl_vector_free ( A->Vold[i][j] );

    for ( i = 0; i < par->var; ++i )
        for ( j = 0; j < par->N; ++j )
            gsl_vector_free ( A->Vnew[i][j] );

    gsl_matrix_free ( A->I );
    gsl_matrix_free ( A->B );

    if ( par->source_flag )
        gsl_matrix_free ( A->Ms );

    printf ( "# SUCCESS # - freed general ADER_DATA.\n" );

}

// ### ### ###
void Free_DG_DATA ( DG_DATA *D, PARA *par ) {
/*
 * D   - DG-Data
 * par - parameter
 *
 * free all allocated arrays and gsl_vectors (sometimes problematic...)
 */

    ASSERT ( D );
    ASSERT ( par );

    int i, k;

// U
    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            gsl_vector_free ( D->U[k][i] );

// pU
    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            gsl_vector_free ( D->pU[k][i] );

// V
    for ( k = 0; k < par->var; ++k )
        free ( D->V[k] );

// beta
    for ( k = 0; k < par->var; ++k )
        free ( D->beta[k] );

// beta_sub
    for ( k = 0; k < par->var; ++k )
        free ( D->beta_sub[k] );

// beta_prev
    for ( k = 0; k < par->var; ++k )
        free ( D->beta_prev[k] );

// Int
    gsl_matrix_free ( D->Int );

// Rec
    gsl_matrix_free ( D->Rec );

// AVG
    free ( D->AVG );

// matrix F
    gsl_matrix_free ( D->F );

// matrix G00
    gsl_matrix_free ( D->G00 );

// matrix G10
    gsl_matrix_free ( D->G10 );

// matrix G01
    gsl_matrix_free ( D->G01 );

// matrix G11
    gsl_matrix_free ( D->G11 );

// matrix H
    if ( par->source_flag )
        gsl_matrix_free ( D->H );

}
