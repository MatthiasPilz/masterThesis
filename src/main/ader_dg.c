#include "ader.h"
#include "main.h"


void ADER_DG_Step ( DG_DATA *D, ADER_DATA *A, PARA *par, PROJ *pro ) {
/*
 * D - contains DG-matrices, old DG-polynomial solution
 * A - contains flux coefficients, DG-predictor
 * par - parameter
 * pro - project related function pointer
 *
 * calculate fully discrete one-step DG-corrector according to eq.(16) in [2]
 * equation is extend with source term! (DG->H matrix)
 */

    ASSERT ( D );
    ASSERT ( A );
    ASSERT ( par );
    ASSERT ( pro );

    int i, k;
    double s;
    gsl_vector *help = gsl_vector_calloc ( par->Mp1 );

    // calc flux and source coefficients
    pro->f_flux_coeff_all ( A, par );

    if ( par->source_flag ) //  pro->f_source_coeff_all ( A, par );
        for ( i = 0; i < par->N; ++i )
            pro->f_source_coeff_i ( A, par, i );


    for ( i = 1; i < par->N-1; ++i ) {  // loop cells

        // calc speed depending on cell average data
        s = fabs ( pro->f_calc_speed ( D->V_grid, par, i ) );

        for ( k = 0; k < par->var; ++k ) {  // loop variables
            gsl_vector_set_zero ( help );   // reset vector

        // G terms f
            gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G11, A->F[k][ i ], 0.0, help );
            gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G10, A->F[k][i+1], 1.0, help );
            gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G00, A->F[k][ i ], 1.0, help );
            gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G01, A->F[k][i-1], 1.0, help );

        // G terms q
            gsl_blas_dgemv ( CblasNoTrans, -s, D->G10, A->Vold[k][i+1], 1.0, help );
            gsl_blas_dgemv ( CblasNoTrans,  s, D->G11, A->Vold[k][ i ], 1.0, help );
            gsl_blas_dgemv ( CblasNoTrans,  s, D->G00, A->Vold[k][ i ], 1.0, help );
            gsl_blas_dgemv ( CblasNoTrans, -s, D->G01, A->Vold[k][i-1], 1.0, help );

            gsl_vector_scale ( help, -0.5*par->CFL );

        // F term
            gsl_blas_dgemv ( CblasNoTrans, par->CFL, D->F, A->F[k][i], 1.0, help );

        // H term (only if source term is needed)
            if ( par->source_flag ) {
                gsl_blas_dgemv ( CblasNoTrans, par->dt, D->H, A->S[k][i], 1.0, help );
            }

            gsl_vector_set_zero ( D->pU[k][i] );                // reset predictor
            gsl_vector_memcpy   ( D->pU[k][i], D->U[k][i] );    // predictor <-- old solution

            gsl_vector_add      ( D->pU[k][i], help );          // predictor += help

        }

    }

    if ( strcmp ( par->boundary, "periodic" ) == 0 ) {
// Periodic Boundaries i = 0
{
    // calc speed depending on cell average data
    s = fabs ( pro->f_calc_speed ( D->V_grid, par, 0 ) );

    for ( k = 0; k < par->var; ++k ) {
        gsl_vector_set_zero ( help );   // reset vector
    // G terms f
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G11, A->F[k][0], 0.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G10, A->F[k][1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G00, A->F[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G01, A->F[k][par->N-1], 1.0, help );

    // G terms q
        gsl_blas_dgemv ( CblasNoTrans, -s, D->G10, A->Vold[k][1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G11, A->Vold[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G00, A->Vold[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -s, D->G01, A->Vold[k][par->N-1], 1.0, help );

        gsl_vector_scale ( help, -0.5*par->CFL );

    // F term
        gsl_blas_dgemv ( CblasNoTrans, par->CFL, D->F, A->F[k][0], 1.0, help );

    // H term
        if ( par->source_flag ) {
            gsl_blas_dgemv ( CblasNoTrans,  par->dt, D->H, A->S[k][0], 1.0, help );
        }

        gsl_vector_set_zero ( D->pU[k][0] );
        gsl_vector_memcpy ( D->pU[k][0], D->U[k][0] );

        gsl_vector_add ( D->pU[k][0], help );
    }
}

// Periodic Boundaries i = par->N-1
{

    // calc speed depending on cell average data
    s = fabs ( pro->f_calc_speed ( D->V_grid, par, par->N-1 ) );

    for ( k = 0; k < par->var; ++k ) {
        gsl_vector_set_zero ( help );   // reset vector
    // G terms f
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G11, A->F[k][par->N-1], 0.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G10, A->F[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G00, A->F[k][par->N-1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G01, A->F[k][par->N-2], 1.0, help );

    // G terms q
        gsl_blas_dgemv ( CblasNoTrans, -s, D->G10, A->Vold[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G11, A->Vold[k][par->N-1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G00, A->Vold[k][par->N-1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -s, D->G01, A->Vold[k][par->N-2], 1.0, help );

        gsl_vector_scale ( help, -0.5*par->CFL );

    // F term
        gsl_blas_dgemv ( CblasNoTrans, par->CFL, D->F, A->F[k][par->N-1], 1.0, help );

    // H term
        if ( par->source_flag ) {
            gsl_blas_dgemv ( CblasNoTrans,  par->dt, D->H, A->S[k][par->N-1], 1.0, help );
        }

        gsl_vector_set_zero ( D->pU[k][par->N-1] );
        gsl_vector_memcpy ( D->pU[k][par->N-1], D->U[k][par->N-1] );

        gsl_vector_add ( D->pU[k][par->N-1], help );
    }
    }
    }

    else if ( strcmp ( par->boundary, "symmetric" ) == 0 || (strcmp ( par->boundary, "hack" ) == 0 )) {
// Periodic Boundaries i = 0
    {

    gsl_vector *temp = gsl_vector_calloc(par->size) ;

    // calc speed depending on cell average data
    s = fabs ( pro->f_calc_speed ( D->V_grid, par, 0 ) );

    for ( k = 0; k < par->var; ++k ) {
        gsl_vector_set_zero ( help );   // reset vector
    // G terms f
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G11, A->F[k][0], 0.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  1.0, D->G10, A->F[k][1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G00, A->F[k][0], 1.0, help );

        int j,t;
        for(t=0;t<par->Mp1;t++)
            for(j=0;j<par->Mp1;j++)
                gsl_vector_set(temp, par->Mp1 - 1 - j + t*par->Mp1, (double)(par->sym_flags[k+par->var]) * gsl_vector_get(A->F[k][0], j + t*par->Mp1));

        gsl_blas_dgemv ( CblasNoTrans, -1.0, D->G01, temp, 1.0, help );

    // G terms q
        gsl_blas_dgemv ( CblasNoTrans, -s, D->G10, A->Vold[k][1], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G11, A->Vold[k][0], 1.0, help );
        gsl_blas_dgemv ( CblasNoTrans,  s, D->G00, A->Vold[k][0], 1.0, help );

        for(t=0;t<par->Mp1;t++)
            for(j=0;j<par->Mp1;j++)
                gsl_vector_set(temp, par->Mp1 - 1 - j + t*par->Mp1, (double)(par->sym_flags[k]) * gsl_vector_get(A->Vold[k][0], j + t*par->Mp1));

        gsl_blas_dgemv ( CblasNoTrans, -s, D->G01, temp, 1.0, help );

        gsl_vector_scale ( help, -0.5*par->CFL );

    // F term
        gsl_blas_dgemv ( CblasNoTrans, par->CFL, D->F, A->F[k][0], 1.0, help );

    // H term
        if ( par->source_flag ) {
            gsl_blas_dgemv ( CblasNoTrans,  par->dt, D->H, A->S[k][0], 1.0, help );
        }

        gsl_vector_set_zero ( D->pU[k][0] );
        gsl_vector_memcpy ( D->pU[k][0], D->U[k][0] );

        gsl_vector_add ( D->pU[k][0], help );
    }

    gsl_vector_free(temp);
    }

// Periodic Boundaries i = par->N-1
    for ( k = 0; k < par->var; ++k )
        gsl_vector_memcpy ( D->pU[k][par->N-1], D->pU[k][par->N-2] );

    }
    else {
        for ( k = 0; k < par->var; ++k ) {
            // i = 0
            gsl_vector_memcpy ( D->pU[k][0], D->pU[k][1] );

            // i = par->N-1
            gsl_vector_memcpy ( D->pU[k][par->N-1], D->pU[k][par->N-2] );
        }
    }

    if ( strcmp ( par->project, "grhdc" ) == 0 )
        grhdc_update_primitives_all(D->pU, par->Mp1, par);

// free temporary vector
    gsl_vector_free ( help );

}

// ### ### ### ### ### ### ###

void Accept_Unlimited_DG ( DG_DATA *DG, FV_DATA *V_AW, WENO_DATA *W_AW, ADER_DATA *A_AW, PARA *par, PARA *par_AW, PROJ *pro ) {
/*
 * DG   - DG DATA
 * V_AW - FV DATA on subcells
 * W_AW - WENO DATA on subcells
 * A_AW - ADER DATA on subcells
 * par  - parameter of grid
 * par_AW - parameter of subgrid
 * pro  - project related function pointer
 *
 * accept predictor solution pU if cell is not troubled
 * otherwise interpolate to subcells and do one ADER WENO step
 */

    ASSERT ( DG );
    ASSERT ( V_AW );
    ASSERT ( W_AW );
    ASSERT ( A_AW );
    ASSERT ( par );
    ASSERT ( par_AW );
    ASSERT ( pro );

    int i, k;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i )      // loop cells
            if ( DG->beta[k][i] != NEED_AW )
                gsl_vector_memcpy ( DG->U[k][i], DG->pU[k][i] );


    // ADER_WENO_subcell
    Interpolation_Prev_Troubled ( DG, V_AW->Avg, DG->Int, par, par_AW );
    Set_Beta_Sub                ( DG->beta, DG->beta_sub, par );
    SUBCELL_WENO_TIME_STEP      ( V_AW, W_AW, A_AW, DG->beta_sub, par_AW, pro );
    Reconstruction              ( V_AW->Avg, DG->U, DG->beta, DG->Rec, par_AW, par );

    if ( strcmp ( par->project, "grhdc" ) == 0 ) {
      grhdc_update_primitives_double_all(V_AW->Avg, par);
      grhdc_update_primitives_all(DG->U, par->Mp1, par);
    }

}






















