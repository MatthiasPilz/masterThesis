#include "ader.h"
#include "ader_weno.h"


void Set_Initial_Iteration_Condition_i ( ADER_DATA *A, gsl_vector ***W, PARA *par, int i ) {
/*
 * A   - Vold are current coefficients for space-time polynomial
 * W   - coefficients for spatial polynomials from WENO reconstruction
 * par - parameter
 * i   - current cell
 *
 * set space-time coefficients as if we have a steady solution --> all time-'layers' have same data in cell i
 */

    ASSERT ( A );
    ASSERT ( W );
    ASSERT ( par );

    int k, j;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( j = 0; j < par->size; ++j )   // loop gaussian quadrature points
            gsl_vector_set ( A->Vold[k][i], j, gsl_vector_get( W[k][i], j%(par->Mp1) ) );

    OUTPUT_EVO ( "set initial iteration conditions for cell" );

}

// ####

double Calc_Change_k ( gsl_vector ***Vold, gsl_vector ***Vnew, int i, int k, PARA *par ) {
/*
 * Vold - space-time coefficients of previous iteration step
 * Vnew - space-time coefficients of current iteration step
 * i    - current cell
 * k    - current variable
 * par  - parameter
 *
 * calculate L2-norm of difference between coefficients of two different iteration steps for a specific cell and variable, return result as double
 */

    ASSERT ( Vold );
    ASSERT ( Vnew );
    ASSERT ( par );
    ASSERT ( i >= 0 && i < par->N );
    ASSERT ( k >= 0 && k < par->var );

    int j;
    double temp = 0.0;

    for ( j = 0; j < par->size; ++ j )  // loop gaussian quadrature points
       temp += pow ( gsl_vector_get(Vold[k][i], j) - gsl_vector_get(Vnew[k][i], j) , 2 );

    // norm accumulated error
    temp = sqrt ( temp );
    temp /= (double)( par->size );

    // return result
    return temp;

}

void Final_Step ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro ) {
        // go to the exact end time

    NEWLINE
    printf ( "*** calculate final small step\n" );
    par->dt = par->tEND - par->time;        printf ( "*** dt  = %.10f\n", par->dt );
                                            printf ( "*** dx  = %.10f\n", par->dx );
    par->CFL = (double)par->dt / par->dx;   printf ( "*** CFL = %.10f\n", par->CFL );

    Set_Beta_Bad            ( beta, par );
    par->iterAvgTemp = par->iterAvg;        // don't add iterAvg from final step.
    par->evalCF = par->evalCountF;
    par->evalPF = par->evalPrediF;
    par->evalCS = par->evalCountS;
    par->evalPS = par->evalPrediS;

    SUBCELL_WENO_TIME_STEP_woStatistics  ( V, W, A, beta, par, pro );

    par->iterAvg = par->iterAvgTemp;
    par->evalCountF = par->evalCF;
    par->evalPrediF = par->evalPF;
    par->evalCountS = par->evalCS;
    par->evalPrediS = par->evalPS;

    par->CFL = par->CFL_sub;
    par->time = par->tEND;

}
