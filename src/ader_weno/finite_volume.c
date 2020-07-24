#include "ader.h"
#include "ader_weno.h"


void Calculate_Flux_Averages ( FV_DATA *V, ADER_DATA *A, PARA *par, int **beta, PROJ *pro ) {
/*
 * V    - contains matrices, flux coefficients and flux cell averages
 * A    - contains current space-time coefficients
 * par  - parameter
 * beta - troubled cell indicator
 * pro  - project related function pointer
 *
 * calculate flux cell averages for all 'troubled' cells according to eq.(4) in [1] using the local Lax-Friedrichs flux eq.(8)
 * matrices are precomputed in mathematica
 */

    ASSERT ( V );
    ASSERT ( A );
    ASSERT ( par );
    ASSERT ( beta );
    ASSERT ( pro );

    int i, k;
    double help;
    double s;

    // calculate f(U) from space-time coefficients only for troubled cells
    for ( i = 0; i < par->N; ++i )    // loop cells TODO ### if ( beta[0][i] != OK ) ###!!!
            pro->f_flux_coeff_i ( A, par, i );


    for ( i = 1; i < par->N; ++i )          // loop cells
        if ( beta[0][i] != OK ) {           // only calculate flux for troubled cells

            s = fabs ( pro->f_calc_speed ( V->Avg, par, i-1 ) );

            for ( k = 0; k < par->var; ++k ) {  // loop variables

                // F(q+,q-) = 0.5 * ( f(q-) + f(q+) - s * ( q+ - q- ) )
                gsl_blas_ddot ( A->F[k][i-1], V->F1_1, &V->Flux[k][i] );
                gsl_blas_ddot ( A->F[k][ i ], V->F1_0, &help );
                    V->Flux[k][i] += help;

                gsl_blas_ddot ( A->Vold[k][i-1], V->F1_1, &help );
                    V->Flux[k][i] += s*help;
                gsl_blas_ddot ( A->Vold[k][ i ], V->F1_0, &help );
                    V->Flux[k][i] -= s*help;

                V->Flux[k][i] *= 0.5;

            }
        }

    // boundary at 0
    if ( strcmp ( par->boundary, "periodic" ) == 0 ) {

        s = fabs ( pro->f_calc_speed ( V->Avg, par, par->N-1 ) );

        for ( k = 0; k < par->var; ++k ) {  // loop variables

            // F(q+,q-) = 0.5 * ( f(q-) + f(q+) - s * ( q+ - q- ) )
            gsl_blas_ddot ( A->F[k][par->N-1], V->F1_1, &V->Flux[k][0] );
            gsl_blas_ddot ( A->F[k][ 0 ], V->F1_0, &help );
                V->Flux[k][0] += help;

            gsl_blas_ddot ( A->Vold[k][par->N-1], V->F1_1, &help );
                V->Flux[k][0] += s*help;
            gsl_blas_ddot ( A->Vold[k][ 0 ], V->F1_0, &help );
                V->Flux[k][0] -= s*help;

            V->Flux[k][0] *= 0.5;

        }
    }
    // boundary at par->xMax
    if ( strcmp ( par->boundary, "periodic" ) == 0 ) {

        s = fabs ( pro->f_calc_speed ( V->Avg, par, par->N-1 ) );

        for ( k = 0; k < par->var; ++k ) {  // loop variables

                gsl_blas_ddot ( A->F[k][par->N-1], V->F1_1, &V->Flux[k][par->N] );
                gsl_blas_ddot ( A->F[k][0], V->F1_0, &help );
                    V->Flux[k][par->N] += help;

                gsl_blas_ddot ( A->Vold[k][par->N-1], V->F1_1, &help );
                    V->Flux[k][par->N] += s*help;
                gsl_blas_ddot ( A->Vold[k][0], V->F1_0, &help );
                    V->Flux[k][par->N] -= s*help;

                V->Flux[k][par->N] *= 0.5;

        }
    }

    // boundary at 0
    if ( strcmp ( par->boundary, "symmetric" ) == 0 || strcmp ( par->boundary, "hack" ) == 0 ) {

        s = fabs ( pro->f_calc_speed ( V->Avg, par, 0 ) );

        gsl_vector *temp = gsl_vector_calloc(par->size);

        for ( k = 0; k < par->var; ++k ) {  // loop variables

            // F(q+,q-) = 0.5 * ( f(q-) + f(q+) - s * ( q+ - q- ) )

            int j,t;
            for(t=0;t<par->Mp1;t++)
                for(j=0;j<par->Mp1;j++)
                    gsl_vector_set(temp, par->Mp1 - 1 - j + t*par->Mp1, (double)(par->sym_flags[k+par->var]) * gsl_vector_get(A->F[k][0], j + t*par->Mp1));

            gsl_blas_ddot ( temp, V->F1_1, &V->Flux[k][0] );
            gsl_blas_ddot ( A->F[k][ 0 ], V->F1_0, &help );
                V->Flux[k][0] += help;

            for(t=0;t<par->Mp1;t++)
                for(j=0;j<par->Mp1;j++)
                    gsl_vector_set(temp, par->Mp1 - 1 - j + t*par->Mp1, (double)(par->sym_flags[k]) * gsl_vector_get(A->Vold[k][0], j + t*par->Mp1));

            gsl_blas_ddot ( temp, V->F1_1, &help );
                V->Flux[k][0] += s*help;
            gsl_blas_ddot ( A->Vold[k][ 0 ], V->F1_0, &help );
                V->Flux[k][0] -= s*help;

            V->Flux[k][0] *= 0.5;

        }

        gsl_vector_free(temp);
    }

}

// ### ### ### ### ### ### ###

void Calculate_Source_Averages ( FV_DATA *V, ADER_DATA *A, PARA *par, int **beta, PROJ *pro ) {
/*
 * V    - contains matrix, source coefficients and source cell averages
 * A    - contains current space-time coefficients
 * par  - parameter
 * beta - troubled cell indicator
 * pro  - project related function pointer
 *
 * calculate source cell averages for all 'troubled' cells according to eq.(7) in [1]
 * matrix are precomputed in mathematica
 */

    ASSERT ( V );
    ASSERT ( A );
    ASSERT ( par );
    ASSERT ( beta );
    ASSERT ( pro );

    int i, k;

    // calculate f(U) from space-time coefficients only for troubled cells
    for ( i = 0; i < par->N; ++i )
        if ( beta[0][i] != OK )
            pro->f_source_coeff_i ( A, par, i );

    for ( i = 0; i < par->N; ++i ) {     // loop cells
        if ( beta[0][i] != OK ) {            // only for troubled cells
            for ( k = 0; k < par->var; ++k ) {   // loop variables
                gsl_blas_ddot ( A->S[k][ i ], V->S1, &V->Source[k][i] );
            }
        }
    }

}

// ### ### ### ### ### ### ###

void Update_Finite_Volume_Scheme ( FV_DATA *V, PARA *par, int **beta ) {
/*
 * V    - contains flux-, source-, and cell averages
 * par  - parameter
 * beta - troubled cell indicator
 *
 * update cell averages using eq.(2) in [1]
 */

    ASSERT ( V );
    ASSERT ( par );
    ASSERT ( beta );

    int k, i;

        if ( par->source_flag ) {   // with source term
            for ( k = 0; k < par->var; ++k ) {  // loop variables
                for ( i = 0; i < par->N; i++ ) {    // loop cells
                    if ( beta[k][i] == NEED_AW ) {      // only troubled cells
                        V->Avg[k][i] = V->Avg[k][i] - par->CFL * ( V->Flux[k][i+1] - V->Flux[k][i] ) + par->dt*V->Source[k][i];
                    }
                }

                    // BOUNDARIES                
                if ( strcmp ( par->boundary, "none" ) == 0 ) {
                    if ( beta[k][0] == NEED_AW )
                    V->Avg[k][0] = V->Avg[k][1];
                    if ( beta[k][par->N-1] == NEED_AW )
                    V->Avg[k][par->N-1] = V->Avg[k][par->N-2];
                }
                else if ( strcmp ( par->boundary, "symmetric" ) == 0 ) {
       // TODO - WHAT IS HAPPENING HERE!?
                    if ( beta[k][par->N-1] == NEED_AW )
                    V->Avg[k][par->N-1] = V->Avg[k][par->N-2];
                }
                else if ( strcmp ( par->boundary, "hack" ) == 0 ) {
                    if ( beta[k][0] == NEED_AW )
                        ;// TODO
                    if ( beta[k][par->N-1] == NEED_AW )
                        ;// TODO
                }

            }
        }
        else {                      // without source term
            for ( k = 0; k < par->var; ++k ) {  // loop variables
                for ( i = 0; i < par->N; i++ ) {    // loop cells
                    if ( beta[k][i] == NEED_AW ) {      // only troubled cells
                        V->Avg[k][i] = V->Avg[k][i] - par->CFL * ( V->Flux[k][i+1] - V->Flux[k][i] );
                    }
                }

                    // BOUNDARIES                
                if ( strcmp ( par->boundary, "none" ) == 0 ) {
                    if ( beta[k][0] == NEED_AW )
                    V->Avg[k][0] = V->Avg[k][1];
                    if ( beta[k][par->N-1] == NEED_AW )
                    V->Avg[k][par->N-1] = V->Avg[k][par->N-2];
                }
                else if ( strcmp ( par->boundary, "symmetric" ) == 0 ) {
       // TODO - WHAT IS HAPPENING HERE!?
                    if ( beta[k][par->N-1] == NEED_AW )
                    V->Avg[k][par->N-1] = V->Avg[k][par->N-2];
                }
                else if ( strcmp ( par->boundary, "hack" ) == 0 ) {
                    if ( beta[k][0] == NEED_AW )
                        ;// TODO
                    if ( beta[k][par->N-1] == NEED_AW )
                        ;// TODO
                }
            }
        }

        if ( strcmp ( par->project, "grhdc" ) == 0 )
            grhdc_update_primitives_double_all(V->Avg,par);

}
