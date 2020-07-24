
#include "ader.h"
#include "main.h"



//// #### #### #### #### ####

void set_initial_q_RK_i ( gsl_vector ****q, gsl_vector ***W, PARA *par, int i ) {
/*
 * q - RK variable
 * W - initial condition from WENO reconstruction
 * par - parameters
 * i - current cell
 *
 * Set q variable with initial condition
 * On each temporal slice values are identical to initial values
 *
 */

    int t, s, k;
    double temp;

    for ( k = 0; k < par->var; ++k ) // variables
        for ( t = 0; t < par->stages; ++t ) // temporal points (stages)
            for ( s = 0; s < par->Mp1; ++s ) {  // spatial points
                temp = gsl_vector_get ( W[k][i], s );
                gsl_vector_set ( q[k][i][t], s, temp );
            }

}

void calc_RK_k ( gsl_vector ****q, gsl_vector ****k, gsl_matrix *D, double *LAM, gsl_vector **F, gsl_vector **S, PARA *par, int i, int t, PROJ *pro ) {
/*
 * q - given values
 * k - store result here
 * par - parameter
 * i - current cell
 * t - current stage
 * pro - project related flux calculation
 *
 * calculate k=g(q) in the t-th step of the RK scheme
 *
 */

    int var;

    pro->f_calc_F_RK ( F, q, par, i, t );

    if ( par->source_flag )
        pro->f_calc_S_RK ( S, q, par, i, t, D, LAM );

    for ( var = 0; var < par->var; ++var ) {
        // flux
        gsl_blas_dgemv ( CblasTrans, -par->CFL, D, F[var], 0.0, k[var][i][t] );

        // source
        if ( par->source_flag ) {
            gsl_vector_scale    ( S[var], par->dt );
            gsl_vector_add      ( k[var][i][t], S[var] );
        }
    }

}


void RK_Predictor ( ADER_DATA *A, gsl_vector ***W, int **beta, PARA *par, PROJ *pro ) {
/*
 * A - will contain predictor
 * W - initial weno reconstruction data
 * par - parameters
 * pro - project related function evaluations
 *
 * in this function the predictor is computed cell-wise using an explizit Runge Kutta scheme for each collocation point in space. (h=1)
 * from the weights k_i we compute the dense output and take the function values at the Gauss Legendre Points as the coefficients of our predictor
 *
 */

    int i, t, s, var, j, k;

    gsl_vector **temp = malloc ( par->var * sizeof ( &temp ) );
    for ( var = 0; var < par->var; ++var )
        temp[var] = gsl_vector_calloc ( par->Mp1 );

    double coeff;

    // calc all k values
    for ( i = 0; i < par->N; ++i ) {    // loop cells
        if ( beta[0][i] != OK ) {           // only for troubled cells

            // initial q for all stages t the same !and all components!
            set_initial_q_RK_i ( A->q_RK, W, par, i );

            // first step
            calc_RK_k ( A->q_RK, A->k_RK, A->D, A->LAM, A->F_RK, A->S_RK, par, i, 0, pro );
            if ( strcmp ( par->project, "grhdc" ) == 0 )
                grhdc_update_primitives_i_t(A->q_RK,par->Mp1,par,i,0);

            // other steps
            for ( t = 1; t < par->stages; ++t ) {
                // calc new q
                for ( var = 0; var < par->var; ++var ) {
                    for ( j = 0; j < t; ++j ) {
                        gsl_vector_memcpy ( temp[var], A->k_RK[var][i][j] );
                        gsl_vector_scale  ( temp[var], A->Aij[t][j] );
                        gsl_vector_add    ( A->q_RK[var][i][t], temp[var] );
                    }
                }

                calc_RK_k ( A->q_RK, A->k_RK, A->D, A->LAM, A->F_RK, A->S_RK, par, i, t, pro );
                if ( strcmp ( par->project, "grhdc" ) == 0 )
                    grhdc_update_primitives_i_t(A->q_RK,par->Mp1,par,i,t);

            }
        }
    }

    // calc dense output --> set A->Vold coefficients
    for ( i = 0; i < par->N; ++i ) {           // loop cells
        for ( var = 0; var < par->var; ++var ) // loop variables
            if ( beta[0][i] != OK ) {           // only for troubled cells
                for ( s = 0; s < par->Mp1; ++s ) {  // loop spatial points
                    for ( t = 0; t < par->Mp1; ++t ) {  // loop time levels
                        coeff = 0.0;

                        for ( k = 0; k < par->stages; ++k )    // loop dense output poly
                            coeff += gsl_vector_get ( A->k_RK[var][i][k],s) * A->Dense[t][k];

                        gsl_vector_set ( A->Vold[var][i], s+(par->Mp1)*t, gsl_vector_get(W[var][i],s) + coeff );
                    }
                }
            }

        if ( beta[0][i] != OK )
            if ( strcmp ( par->project, "grhdc" ) == 0 )
                grhdc_update_primitives_i(A->Vold, par->size, par, i);


    }

  /// FREE TEMPORARY VECTOR
    for ( var = 0; var < par->var; ++var )
        gsl_vector_free ( temp[var] );

    free ( temp );

}













