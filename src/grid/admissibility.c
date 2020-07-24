#include "ader.h"
#include "grid.h"



static double voronoi_neigh_MIN ( double *V, int cur, PARA *par ) {
/*
 * find minimum cell average in voronoi neighbouring cells
 */

    ASSERT ( cur >= 0 && cur < par->N );

    double min = INFINITY;
    int count;

    for ( count = (cur-1)*par->NSC; count < (cur+2)*par->NSC; ++count )
        if ( V[count] < min )
            min = V[count];

    return min;

}

static double voronoi_neigh_MAX ( double *V, int cur, PARA *par ) {
/*
 * find maximum cell average in voronoi neighbouring cells
 */

    ASSERT ( cur >= 0 && cur < par->N );

    double max = -INFINITY;
    int count;

    for ( count = (cur-1)*par->NSC; count < (cur+2)*par->NSC; ++count )
        if ( V[count] > max )
            max = V[count];

    return max;

}

void Numerical_Admissibility ( DG_DATA *DG, double **V, PARA *par ) {
/*
 * DG  - DG_DATA
 * V   - cell averages on subcells
 * par - parameter on subcells
 *
 * calculate numerical admissibility according to eq.(25) in [2]
 */

    ASSERT ( DG );
    ASSERT ( V );
    ASSERT ( par );

    int i, k;
    int cur;
    double min, max;
    double delta;
    double delta0 = 1e-4;

    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = par->NSC; i < par->N-par->NSC; ++i ) { // loop subcells but not boundary
            // find voronoi cell:
            cur = (int)i / ((int)par->NSC);

            if ( DG->beta[k][cur] != NEED_AW ) {
                // calculate minimum and maximum value of variable
                // in voronoi neighbouring cells
                min = voronoi_neigh_MIN ( V[k], cur, par );
                max = voronoi_neigh_MAX ( V[k], cur, par );

                delta = fabs ( 1e-3 * ( max - min ) );

/// modified version (see Zanotti et al. '15)
                if ( delta < delta0 )
                    delta = delta0;

                // set troubled cell indicator
                if ( (min - delta) > DG->V[k][i] ) {
                    DG->beta[k][ cur ] = NEED_AW;
                    if ( DG->beta[k][cur-1] != NEED_AW )
                        DG->beta[k][cur-1] = AW_NEIGH;

                    if ( DG->beta[k][cur+1] != NEED_AW )
                        DG->beta[k][cur+1] = AW_NEIGH;
                }
                else if ( (max + delta ) < DG->V[k][i] ) {
                    DG->beta[k][ cur ] = NEED_AW;
                    if ( DG->beta[k][cur-1] != NEED_AW )
                        DG->beta[k][cur-1] = AW_NEIGH;

                    if ( DG->beta[k][cur+1] != NEED_AW )
                        DG->beta[k][cur+1] = AW_NEIGH;
                }
            }

        }


}

// ### ### ### ### ### ### ###

static void calc_boundary_polynomial_values ( DG_DATA *D, double *Um, double *Up, PARA *par ) {
// calculate upper and lower boundary value of predictor polynomial

    int i, k, j;
    double tempLow = 0.0;
    double tempUp  = 0.0;

    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i ) {

            tempLow = 0.0;
            tempUp  = 0.0;

            for ( j = 0; j < par->Mp1; ++j ) {
                tempLow += D->BLow[j] * gsl_vector_get ( D->pU[k][i], j );
                tempUp  += D->BUp[j]  * gsl_vector_get ( D->pU[k][i], j );
            }

            Um[i] = tempLow;
            Up[i] = tempUp;
        }

}

static double MinMod ( double a1, double a2, double a3 ) {
// Qiu Shu Eq. (2.4)

    double temp = 1e10;

    if ( (SIGNOF (a1) == SIGNOF (a2)) && (SIGNOF (a2) == SIGNOF (a3)) ) {
        if ( fabs(a1) < temp ) temp = fabs(a1);
        if ( fabs(a2) < temp ) temp = fabs(a2);
        if ( fabs(a3) < temp ) temp = fabs(a3);

        return ( SIGNOF(a1) * temp );
    }
    else
        return 0.0;

}

static double ModMinMod ( double a1, double a2, double a3, double MH2 ) {
// Qiu Shu Eq. (2.5)

    if ( fabs(a1) <= MH2 )
        return a1;
    else
        return ( MinMod ( a1, a2, a3 ) );

}

void Troubled_Cell_Indicator ( DG_DATA *DG, PARA *par ) {
/*
 * DG - contains predictor pU and cell averages V which is needed for troubled cell indicator
 * par - parameter
 *
 * find troubled cells according to algorithm of Cockburn & Shu!? applying minmod function
 */

    int i, k;
    double ut, utt, dpu, dmu; //u-tilde, u-tilde-tilde, delta+u, delta-u
    double utmod, uttmod;       //u-tilde-mod, u-tilde-tilde-mod

    // not yet optimized...
    double *Um = calloc ( par->N, sizeof ( double ) );  // lower boundary values
    double *Up = calloc ( par->N, sizeof ( double ) );  // upper boundary values

    calc_boundary_polynomial_values ( DG, Um, Up, par );

    for ( k = 0; k < par->var; ++k ) {
        for ( i = 1; i < par->N-1; ++i ) {

            // calc all variables
            ut  = Um[i+1] - DG->V_grid_pre[k][i];
            utt = DG->V_grid_pre[k][i] - Up[i-1];
            dpu = DG->V_grid_pre[k][i+1] - DG->V_grid_pre[k][i];
            dmu = DG->V_grid_pre[k][i] - DG->V_grid_pre[k][i-1];

            // find minmod variables
        /// standard minmod
//            utmod  = MinMod ( ut,  dpu, dmu );
//            uttmod = MinMod ( utt, dpu, dmu );

        /// modified minmod
            utmod  = ModMinMod ( ut,  dpu, dmu, par->MH2 );
            uttmod = ModMinMod ( utt, dpu, dmu, par->MH2 );

            // decide whether cell is troubled or not and set neighbours to AW_NEIGH
            if ( utmod != ut || uttmod != utt ) {
                DG->beta[k][i] = NEED_AW;

                if ( DG->beta[k][i-1] != NEED_AW )
                    DG->beta[k][i-1] = AW_NEIGH;

                if ( DG->beta[k][i+1] != NEED_AW )
                    DG->beta[k][i+1] = AW_NEIGH;
            }

        }

    }

    free ( Um );
    free ( Up );

}
















