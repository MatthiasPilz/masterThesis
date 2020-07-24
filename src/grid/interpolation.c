#include "ader.h"
#include "grid.h"


void Interpolation ( gsl_vector ***U, double **V, gsl_matrix *Int, PARA *parU, PARA *parV ) {
/*
 * U    - DG-polynomial
 * V    - going to contain cell averages
 * Int  - precomputed interpolation matrix
 * parU - parameter on grid
 * parV - parameter on subgrid
 *
 * calculate interpolation of DG-polynomial on subcells according to eq.(18) in [2]
 */

    ASSERT ( U );
    ASSERT ( V );
    ASSERT ( parU );
    ASSERT ( parV );

    // temporary vector
    gsl_vector *help = gsl_vector_calloc ( parU->NSC );

    int i, k, j;
    for ( k = 0; k < parU->var; ++k )   // loop variables
        for ( i = 0; i < parU->N; ++ i ) {  // loop grid cells
            // matrix multiplication
            gsl_blas_dgemv ( CblasNoTrans, 1.0, Int, U[k][i], 0.0, help );

            // store data in subcell array
            for ( j = 0; j < parU->NSC; ++j )   // loop subcells in cell
                V[k][i*parU->NSC + j] = gsl_vector_get ( help, j );

        }

   gsl_vector_free ( help );

}

// #### #### ####

void Interpolation_Prev_Troubled ( DG_DATA *DG, double **V, gsl_matrix *Int, PARA *parU, PARA *parV ) {

    ASSERT ( DG );
    ASSERT ( V );
    ASSERT ( parU );
    ASSERT ( parV );

    gsl_vector *help = gsl_vector_calloc ( parU->NSC );

    int i, k, j;
    for ( k = 0; k < parU->var; ++k )
        for ( i = 0; i < parU->N; ++ i ) {
            if ( (DG->beta_prev[k][i] == OK) && (DG->beta[k][i] == NEED_AW) ) {

                gsl_blas_dgemv ( CblasNoTrans, 1.0, Int, DG->U[k][i], 0.0, help );

                for ( j = 0; j < parU->NSC; ++j )
                    V[k][i*parU->NSC + j] = gsl_vector_get ( help, j );
            }
        }

   gsl_vector_free ( help );

}

// #### #### ####

void Interpolation_Good_Cells ( DG_DATA *DG, double **V, gsl_matrix *Int, PARA *parU, PARA *parV ) {

    ASSERT ( DG );
    ASSERT ( V );
    ASSERT ( parU );
    ASSERT ( parV );

    gsl_vector *help = gsl_vector_calloc ( parU->NSC );

    int i, k, j;
    for ( k = 0; k < parU->var; ++k )
        for ( i = 0; i < parU->N; ++ i ) {
            if ( (DG->beta_prev[k][i] != NEED_AW) && (DG->beta[k][i] != NEED_AW) ) {

                gsl_blas_dgemv ( CblasNoTrans, 1.0, Int, DG->U[k][i], 0.0, help );

                for ( j = 0; j < parU->NSC; ++j )
                    V[k][i*parU->NSC + j] = gsl_vector_get ( help, j );
            }
        }

   gsl_vector_free ( help );

}

// #### #### ####

void Reconstruction ( double **V, gsl_vector ***U, int **beta, gsl_matrix *Rec, PARA *parV, PARA *parU ) {

    ASSERT ( U );
    ASSERT ( V );
    ASSERT ( parU );
    ASSERT ( parV );

    gsl_vector *help = gsl_vector_calloc ( parU->NSC );

    int i, k, j;

    for ( k = 0; k < parV->var; ++k )
        for ( i = 0; i < parU->N; ++ i ) {
            if ( beta[k][i] == NEED_AW ) {

                for ( j = 0; j < parU->NSC; ++j )
                    gsl_vector_set ( help, j, V[k][i*parU->NSC + j] );

                gsl_blas_dgemv ( CblasNoTrans, 1.0, Rec, help, 0.0, U[k][i] );

            }
        }

   gsl_vector_free ( help );

}

void Cell_Averages ( gsl_vector ***U, double **V, double *AVG, PARA *par ) {

    int i, k, j;

    for ( k = 0; k < par->var; ++k ) {
        for ( i = 0; i < par->N; ++i ) {
            V[k][i] = 0.0;

            for ( j = 0; j < par->Mp1; ++j ) {
                V[k][i] += AVG[j] * gsl_vector_get (U[k][i], j);
            }
        }
    }

}
