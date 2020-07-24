#include "ader.h"
#include "data.h"

void Init_FV_DATA_general ( FV_DATA *V, PARA *par ) {
/*
 * V   - FV_DATA
 * par - parameter
 *
 * initialize precomputed flux and source vectors
 * (todo: use function for reading vectors similar to matrices...)
 */

    int i = 0;
    double temp;
    int check;
    char filename[255];

// first flux vector
    sprintf ( filename, "../math/%d/F1_0", par->M );

    FILE *data1 = fopen ( filename, "r" );	// open file for reading
    ASSERT ( data1 );

    while ( !feof( data1 ) ) {			// go through file number by number
        check = fscanf ( data1, "%lf", &temp);
        ASSERT ( check );
        gsl_vector_set ( V->F1_0, i++, temp );
    }
    fclose ( data1 );


// second flux vector
    i = 0;
    sprintf ( filename, "../math/%d/F1_1", par->M );

    FILE *data2 = fopen ( filename, "r" );	// open file for reading
    ASSERT ( data2 );

    while ( !feof( data2 ) ) {			// go through file number by number
        check = fscanf ( data2, "%lf", &temp);
        ASSERT ( check );
        gsl_vector_set ( V->F1_1, i++, temp );
    }
    fclose ( data2 );


    if ( par->source_flag ) {
// source vector
        i = 0;
        sprintf ( filename, "../math/%d/S1", par->M );

        FILE *data = fopen ( filename, "r" );	// open file for reading
        ASSERT ( data );

        while ( !feof( data ) ) {			// go through file number by number
            check = fscanf ( data, "%lf", &temp);
            ASSERT ( check );
            gsl_vector_set ( V->S1, i++, temp );
        }
        fclose ( data );
    }

    printf ( "# SUCCESS # - initialized FV_DATA struct.\n" );

}

// ### ### ### ### ### ### ###

static void init_rec_matrices ( gsl_matrix **Rec, PARA *par ) {

    ASSERT ( Rec );
    ASSERT ( par );
    ASSERT ( par->n_stencil > 0 );

    int i;
    char fullPath[255];

    for ( i = 0; i < par->n_stencil; ++i ) {
        sprintf ( fullPath, "../math/%d/WENO_Rec%d", par->M, i );
        Read_Gsl_Matrix ( Rec[i], par->Mp1, par->Mp1, fullPath );
    }

    printf ( "# SUCCESS # - initialized all reconstruction matrices.\n" );

}

static void init_lambda ( double *L, PARA *par ) {

    ASSERT ( L );
    ASSERT ( par );
    ASSERT ( par->n_stencil > 0 );

    if ( par->n_stencil == 2 )  {
        L[0] = 1.0;
        L[1] = par->WENO_lam;
    }

    if ( par->n_stencil == 3 ) {
///         ideal weights
//        L[0] = 0.3;
//        L[1] = 0.6;
//        L[2] = 0.1;

///         Dumbser's weights
        L[0] = 1.0;
        L[1] = par->WENO_lam;
        L[2] = 1.0;
    }
    else if ( par->n_stencil == 4 ) {
        L[0] = 1.0;
        L[1] = par->WENO_lam;
        L[2] = par->WENO_lam;
        L[3] = 1.0;
    }
    else
        printf ( "*** ERROR - init_lambda - n_stencil wrong\n" );

    printf ( "# SUCCESS # - initialized lambda.\n" );

}

static void init_OIM ( double **OIM, PARA *par ) {

    ASSERT ( OIM );
    ASSERT ( par );

    char fullPath[255];
    sprintf ( fullPath, "../math/%d/WENO_OIM", par->M );

    Read_Double_Matrix      ( OIM, par->Mp1, par->Mp1, fullPath );

    printf ( "# SUCCESS # - initialized OIM matrix.\n" );

}

static void init_shift ( int *shift, PARA *par ) {

    int i;

    shift[0] = 0;
    shift[par->n_stencil-1] = par->M;

    if ( par->n_stencil == 3 ) {
        shift[1] = (int)(par->M / 2);
    }
    else if ( par->n_stencil == 4 ) {
        shift[1] = (int)(par->M / 2);
        shift[2] = (int)(par->M / 2) + 1;
    }
    else
        printf ( "*** ERROR - init_shift - n_stencil wrong\n" );

}

void Init_WENO_DATA_general ( WENO_DATA *W, PARA *par ) {
/*
 * W   - WENO_DATA
 * par - parameter
 *
 * initialize precomputed matrices and shift vector
 */

    ASSERT ( W );
    ASSERT ( par );

    init_rec_matrices ( W->Rec, par );
    init_lambda ( W->l_s, par );
    init_OIM ( W->OIM, par );
    init_shift ( W->shift, par );

    printf ( "# SUCCESS # - initialized WENO_DATA struct.\n" );

}

// ### ### ### ### ### ### ###

void Init_ADER_DATA_general ( ADER_DATA *A, PARA *par ) {
/*
 * A   - ADER_DATA
 * par - parameter
 *
 * initialize precomputed matrices
 */

    ASSERT ( A );
    ASSERT ( par );

    char fullPath[255];

    sprintf ( fullPath, "../math/%d/ADER_I", par->M );
    Read_Gsl_Matrix ( A->I, par->size, par->Mp1, fullPath );

    sprintf ( fullPath, "../math/%d/ADER_B", par->M );
    Read_Gsl_Matrix ( A->B, par->size, par->size, fullPath );

    if ( par->source_flag ) {
        sprintf ( fullPath, "../math/%d/ADER_M", par->M );
        Read_Gsl_Matrix ( A->Ms, par->size, par->size, fullPath );
    }

    // LAM
    sprintf ( fullPath, "../math/%d/LAM", par->M );
    Read_Double_Vector ( A->LAM, par->Mp1, fullPath );

    // matrix D
    sprintf ( fullPath, "../math/%d/DerivLagr", par->M );
    Read_Gsl_Matrix ( A->D, par->Mp1, par->Mp1, fullPath );

    if ( strcmp ( par->scheme, "RK" ) == 0 ) {

        // Butcher Aij Matrix
        sprintf ( fullPath, "../math/%d/RK_%d_Aij", par->M, par->stages );
        Read_Double_Matrix ( A->Aij, par->stages, par->stages, fullPath );

        // Dense Output
        sprintf ( fullPath, "../math/%d/Dense_Output%d", par->M, par->stages );
        Read_Double_Matrix ( A->Dense, par->Mp1, par->stages, fullPath );
    }

    printf ( "# SUCCESS # - initialized ADER_DATA struct.\n" );

}

// ### ### ### ### ### ### ###


void Init_DG_DATA_general ( FV_DATA *V, WENO_DATA *W, DG_DATA *DG, PARA *par, PROJ *pro ) {
/*
 * V   - FV_DATA
 * W   - WENO_DATA
 * DG  - DG_DATA
 * par - parameter
 * pro - project related functions
 *
 * construct DG-polynomials from WENO reconstruction of given cell averages
 * initialize precomputed matrices
 */

    ASSERT ( V );
    ASSERT ( W );
    ASSERT ( DG );
    ASSERT ( par );

    int i;

    Set_Beta_Bad ( DG->beta, par );
    Set_Beta_Bad ( DG->beta_prev, par );

    // WENO reconstruction
    FV_2_WENO_DATA          ( V, W, par, DG->beta );
    Solve_Algebraic_System  ( W, par, DG->beta );
    Calculate_WENO_Coeff    ( W, par, DG->beta );

    // store WENO coefficients as coefficients for DG-polynomial
    int j, k;
    for ( k = 0; k < par->var; ++k )
        for ( i = 0; i < par->N; ++i )
            for ( j = 0; j < par->Mp1; ++j )
                gsl_vector_set ( DG->U[k][i], j, gsl_vector_get ( W->WT[k][i], j ) );

    // Int
    char fullPath[1023];
    sprintf ( fullPath, "../math/%d/Interpol", par->M );
    Read_Gsl_Matrix ( DG->Int, par->NSC, par->Mp1, fullPath );

    // Rec
    sprintf ( fullPath, "../math/%d/Reconstr", par->M );
    Read_Gsl_Matrix ( DG->Rec, par->Mp1, par->NSC, fullPath );

    // AVG
    sprintf ( fullPath, "../math/%d/AVG", par->M );
    Read_Double_Vector ( DG->AVG, par->Mp1, fullPath );

    // LAM
    sprintf ( fullPath, "../math/%d/LAM", par->M );
    Read_Double_Vector ( DG->LAM, par->Mp1, fullPath );

    // BLow
    sprintf ( fullPath, "../math/%d/BoundLow", par->M );
    Read_Double_Vector ( DG->BLow, par->Mp1, fullPath );

    // Bup
    sprintf ( fullPath, "../math/%d/BoundHigh", par->M );
    Read_Double_Vector ( DG->BUp, par->Mp1, fullPath );


    // DG matrices
    sprintf ( fullPath, "../math/%d/DG_F", par->M );
    Read_Gsl_Matrix ( DG->F, par->Mp1, par->size, fullPath );

    sprintf ( fullPath, "../math/%d/DG_G00", par->M );
    Read_Gsl_Matrix ( DG->G00, par->Mp1, par->size, fullPath );

    sprintf ( fullPath, "../math/%d/DG_G01", par->M );
    Read_Gsl_Matrix ( DG->G01, par->Mp1, par->size, fullPath );

    sprintf ( fullPath, "../math/%d/DG_G10", par->M );
    Read_Gsl_Matrix ( DG->G10, par->Mp1, par->size, fullPath );

    sprintf ( fullPath, "../math/%d/DG_G11", par->M );
    Read_Gsl_Matrix ( DG->G11, par->Mp1, par->size, fullPath );

    if ( par->source_flag ) {
        sprintf ( fullPath, "../math/%d/DG_H", par->M );
        Read_Gsl_Matrix ( DG->H, par->Mp1, par->size, fullPath );
    }

    Set_Beta_Bad ( DG->beta, par );
    Set_Beta_Sub ( DG->beta, DG->beta_sub, par );

    printf ( "# SUCCESS # - initialized DG_DATA struct.\n" );


}




























