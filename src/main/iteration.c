#include "ader.h"
#include "main.h"


void Iteration ( ADER_DATA *A, gsl_vector ***W, int **beta, PARA *par, PROJ *pro ) {
/*
 * A    - ADER DATA
 * W    - coefficients of polynomial
 * par  - parameter
 * beta - troubled cell indicator
 * pro  - project related function pointer
 *
 * calculate local space-time predictor from given polynomial in cell - see eq.(41) in [1]
 * save final space-time coefficients in A->Vold
 */

    ASSERT ( A );
    ASSERT ( W );
    ASSERT ( par );
    ASSERT ( beta );
    ASSERT ( pro );

    int i, k;
    int iter = 0;
    double change;
    double change_k;
    double average_iter = 0.0;
    double average_chan = 0.0;

    for ( i = 0; i < par->N; ++i ) {    // loop cells
        iter = 0;
        change = 10000;

        if ( beta[0][i] != OK ) {           // only for troubled cells

            Set_Initial_Iteration_Condition_i ( A, W, par, i );

            do {    // iteration loop
                pro->f_flux_coeff_i ( A, par, i );
                if ( par->source_flag )
                    pro->f_source_coeff_i ( A, par, i );

                change_k = 0.0;
                for ( k = 0; k < par->var; ++k ) {  // loop variables

                    gsl_blas_dgemv ( CblasNoTrans, -par->CFL, A->B, A->F[k][i], 0.0, A->Vnew[k][i] );
                    gsl_blas_dgemv ( CblasNoTrans, 1.0, A->I, W[k][i], 1.0, A->Vnew[k][i] );
                    if ( par->source_flag ) {
                        gsl_blas_dgemv ( CblasNoTrans, par->dt, A->Ms, A->S[k][i], 1.0, A->Vnew[k][i] );
                    }

                }

                if ( strcmp ( par->project, "grhdc" ) == 0 )
                    grhdc_update_primitives_i(A->Vnew,par->size,par,i);

                for ( k = 0; k < par->var; ++k ) {  // loop variables
                    change_k += Calc_Change_k ( A->Vold, A->Vnew, i, k, par );
                    gsl_vector_swap ( A->Vold[k][i], A->Vnew[k][i] );
                }

                change = change_k;
                iter++;
            } while ( iter < par->iterMAX && change > par->eps );

            // save iteration statistic... todo: switch on/off via flag!?
            average_iter += iter;
            average_chan += change;

//            if ( change > par->changeMAX ) {par->changeMAX = change;}
//            if ( par->iterFail && (change >= 10*par->eps) ) {par->iterFail = 1;}
        }

    }

    // statistic...
    par->iterAvg += ( average_iter / (1.0*par->N) );
    par->chanAvg += ( average_chan / (1.0*par->N) );

}
