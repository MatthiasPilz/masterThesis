#include "ader.h"
#include "ader_weno.h"

void Solve_Algebraic_System ( WENO_DATA *W, PARA *par, int **beta ) {
/*
 * W - WENO Data struct contains Matrices and current data
 * par - parameter
 *
 * Calculate WENO coefficients according to eq.(14)
 */

    ASSERT ( W );
    ASSERT ( par );

    int k, i, s;

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i ) {

            if ( beta[k][i] != OK ) {
                for ( s = 0; s < par->n_stencil; ++s )
                    gsl_blas_dgemv( CblasNoTrans, 1.0, W->Rec[s], W->gslData[k][i][s], 0.0, W->gslWs[k][i][s] );
            }
        }

}

// ### ### ### ### ### ### ###

static void calc_sigma ( double **sigma, double **OIM, gsl_vector ***w, PARA *par, int *beta ) {
    // sigma - empty matrix which is going to store the oscillation indicator for each cell and stencil
    // OIM - precomputed oscillation indicator matrix according to eq.(19)
    // w - solution to algebraic system eq.(15)
    // M - order of scheme
    // N - number of cells on interval [0,1]
// calculate sigma according to eq.(18)

    ASSERT ( sigma );
    ASSERT ( OIM );
    ASSERT ( w );
    ASSERT ( par->N > 2*par->Mp1 );

    int i, s, p, m;
    double sum;

    for ( i = 0; i < par->N; i++ ) {            // loop over cells
        if ( beta[i] != OK ) {
            for ( s = 0; s < par->Mp1; s++ ) {			// loop over stencils
                sum = 0.0;

                for ( p = 0; p < par->Mp1; p++ ) 			// sum over p and m
                    for ( m = 0; m < par->Mp1; m++ )
                        sum += ( OIM[p][m] * gsl_vector_get(w[i][s], p) * gsl_vector_get(w[i][s], m) );

                sigma[i][s] = sum;
//            printf ( "i = %d\ts = %d\tsigma = %.3e\n", i, s, sigma[i][s] );
            }
        }
    }
    OUTPUT_EVO ( "calculated sigma." );
}

static void calc_w_tilde ( double **w_tilde, double *l_s, double **sigma, PARA *par, int *beta ) {
    // w_tilde - is going to contain solution
    // l_s - lambda coefficients
    // sigma - oscillation indicator
    // M - order of scheme
    // N - number of cells on interval [0,1]
// calculate w_tilde according to eq.(17)

    ASSERT ( w_tilde );
    ASSERT ( l_s );
    ASSERT ( sigma );
    ASSERT ( par->N > 2*par->M+1 );

    int i, s;

    int r = par->WENO_R;		// given value for exponent was 8 !!! unstable!?
    double eps = par->WENO_eps; // given value for epsilon summand was 1e-14

    for ( i = 0; i < par->N; i++ ) 			// loop over cells
        if ( beta[i] != OK ) {
            for ( s = 0; s < par->Mp1; s++ )			// loop over stencils
                w_tilde[i][s] = l_s[s] / ( pow ( ( sigma[i][s] + eps ), r ) );
        }

    OUTPUT_EVO ( "calculated w_tilde." );
}

static void calc_w_s ( double **w_s, double **w_tilde, PARA *par, int *beta ) {
    // w_s - is going to contain solution
    // w_tilde - special weight according to eq.(17)
    // M - order of scheme
    // N - number of cells on interval [0,1]
// calculate w_s according to eq.(17)

    ASSERT ( w_s );
    ASSERT ( w_tilde );
    ASSERT ( par->N > 2*par->M+1 );

    int i, s, k;
    double temp;

    for ( i = 0; i < par->N; i++ ) {			// loop over cells
        if ( beta[i] != OK ) {
            temp = 0.0;
            for ( k = 0; k < par->Mp1; k++ )			// calc sum of weights
                temp += w_tilde[i][k];

            for ( s = 0; s < par->Mp1; s++ )	{	// eq.(17) for each stencil
                w_s[i][s] = w_tilde[i][s] / temp;
//            printf ( "i = %d\ts = %d\tws = %.10e\n", i, s, w_s[i][s] );
            }
//        NEWLINE
        }
    }

    OUTPUT_EVO ( "calculated w_s." );
}

static void calc_coeff ( gsl_vector **final, double **w_s, gsl_vector ***w, PARA *par, int *beta ) {
    // final - is going to contain WENO-coefficients
    // w_s - weights according to eq.(17)
    // w - solution to algebraic system eq.(15)
    // M - order of scheme
    // N - number of cells on interval [0,1]
// calculate WENO coefficients according to eq.(16)

    ASSERT ( final );
    ASSERT ( w_s );
    ASSERT ( w );
    ASSERT ( par->N > 2*par->M+1 );

    int i, s, p;
    double temp;

    for ( i = 0; i < par->N; i++ )			// loop over cells
        if ( beta[i] != OK ) {
            for ( p = 0; p < par->Mp1; p++ ) {			// loop over points in stencil
                temp = 0.0;

                for ( s = 0; s < par->Mp1; s++ )			// sum over stencil
                    temp += ( w_s[i][s] * gsl_vector_get( w[i][s], p ) );

                gsl_vector_set ( final[i], p, temp );
            }
        }

    OUTPUT_EVO ( "calculated WENO coefficients." );
}

void Calculate_WENO_Coeff ( WENO_DATA *W, PARA *par, int **beta ) {

    ASSERT ( W );
    ASSERT ( par );

    int k;
    for ( k = 0; k < par->var; ++k ) {
        calc_sigma      ( W->sigma, W->OIM, W->gslWs[k], par, beta[k] );
        calc_w_tilde    ( W->w_tilde, W->l_s, W->sigma, par, beta[k] );
        calc_w_s        ( W->w_s, W->w_tilde, par, beta[k] );
        calc_coeff      ( W->WT[k], W->w_s, W->gslWs[k], par, beta[k] );
    }

}

// ### ### ### ### ### ### ###




