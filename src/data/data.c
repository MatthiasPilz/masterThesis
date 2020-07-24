#include "ader.h"
#include "data.h"

void Read_Gsl_Matrix ( gsl_matrix *A, const int sizeR, const int sizeC, char *fullPath ) {
/*
 * A - resulting matrix
 * sizeR - number of rows
 * sizeC - number of columns
 * fullPath - path from .../exe/ to output file of mathematica
 *
 * read matrix from file (mathematica output) and store as gsl_matrix
 */

    ASSERT ( A );
    ASSERT ( sizeR > 0 );
    ASSERT ( sizeC > 0 );
    ASSERT ( fullPath );

    int i = 0;
    int column, row;
    double temp;
    int check;

    FILE *fp = fopen ( fullPath, "r" ); // open file
    ASSERT ( fp );

    while ( !feof ( fp ) ) {            // go through file till the end
        check = fscanf ( fp, "%lf", &temp);     // read line
        ASSERT ( check );

        // calculate 2d-index, i -> (row,column)
        column = i % sizeC;
        row = (int) ( i - column ) / sizeC;

        ASSERT ( column >= 0 && column < sizeC );
        ASSERT ( row >= 0 && row < sizeR );

        gsl_matrix_set ( A, row, column, temp );
        i++;
    }
    fclose ( fp );

    printf ( "# SUCCESS # - read gsl-matrix: %s\n", fullPath );

}

// ### ### ### ### ### ### ###

void Read_Double_Matrix ( double **A, const int sizeR, const int sizeC, char *fullPath ) {
/*
 * A - resulting matrix
 * sizeR - number of rows
 * sizeC - number of columns
 * fullPath - path from .../exe/ to output file of mathematica
 *
 * read matrix from file (mathematica output) and store as double array
 */

    ASSERT ( A );
    ASSERT ( sizeR > 0 );
    ASSERT ( sizeC > 0 );
    ASSERT ( fullPath );

    int i = 0;
    int column, row;
    double temp;
    int check;

    FILE *fp = fopen ( fullPath, "r" );
    ASSERT ( fp );

    while ( !feof ( fp ) ) {
        check = fscanf ( fp, "%lf", &temp);
        ASSERT ( check );

        column = i % sizeC;
        row = (int) ( i - column ) / sizeC;

        ASSERT ( column >= 0 && column < sizeC );
        ASSERT ( row >= 0 && row < sizeR );

        A[row][column] = temp;
        i++;
    }
    fclose ( fp );

    printf ( "# SUCCESS # - read gsl-matrix: %s\n", fullPath );

}

// ### ### ### ### ### ### ###

void Read_Double_Vector ( double *A, const int size, char *fullPath ) {
/*
 * A - resulting vector
 * size - number of entries
 * fullPath - path from ../exe/ to output file of mathematica
 *
 * read vector from file and store as double array
 */

    int i = 0;
    double temp;
    double check;

    FILE *data = fopen ( fullPath, "r" );   // open file for reading
    ASSERT ( data );

    while ( !feof( data ) ) {	// go through file line by line
        check = fscanf ( data, "%lf", &temp);
        ASSERT ( check );
        A[i++] = temp;		// write result in array and increment i
    }
    fclose ( data );			// close file

    ASSERT ( i == size );
}

// ### ### ### ### ### ### ###

void Read_GSL_Vector ( gsl_vector *A, const int size, char *fullPath ) {
/*
 * A - resulting gsl vector
 * size - number of entries
 * fullPath - path from ../exe/ to output file of mathematica
 *
 * read vector from file and store as gsl vector
 */

    int i = 0;
    double temp;
    double check;

    FILE *data = fopen ( fullPath, "r" );   // open file for reading
    ASSERT ( data );

    while ( !feof( data ) ) {	// go through file line by line
        check = fscanf ( data, "%lf", &temp);
        ASSERT ( check );
        gsl_vector_set ( A, i, temp );
        i++;
    }
    fclose ( data );			// close file

    ASSERT ( i == size );
}

// ### ### ### ### ### ### ###

static int translate_boundary ( int boundary, int N ) {
/*
 * boundary - current cell that is not on 'real' grid
 * N - number of 'real' cells
 *
 * use periodic boundary conditions to wrap boundary back on 'real' grid
 * return index
 *
 * +- 10 should be more than sufficient. todo: generalise to arbitrary orders...
 */

    int result;

    if ( boundary ==  N  )	result = 0;
    else if ( boundary == -1 ) 	result = N-1;
    else if ( boundary == N+1 )	result = 1;
    else if ( boundary == -2 ) 	result = N-2;
    else if ( boundary == N+2 )	result = 2;
    else if ( boundary == -3 ) 	result = N-3;
    else if ( boundary == N+3 )	result = 3;
    else if ( boundary == -4 ) 	result = N-4;
    else if ( boundary == N+4 )	result = 4;
    else if ( boundary == -5 )  result = N-5;
    else if ( boundary == N+5 )	result = 5;
    else if ( boundary == -6 )  result = N-6;
    else if ( boundary == N+6 )	result = 6;
    else if ( boundary == -7 )  result = N-7;
    else if ( boundary == N+7 )	result = 7;
    else if ( boundary == -8 )  result = N-8;
    else if ( boundary == N+8 )	result = 8;
    else if ( boundary == -9 )  result = N-9;
    else if ( boundary == N+9 )	result = 9;
    else if ( boundary == -10 )  result = N-10;
    else if ( boundary == N+10 ) result = 10;
    else {
        printf ( "*** ERROR in D2GV ( %d )\n", boundary );
        result = -1;
    }

    return result;

}

void FV_2_WENO_DATA ( FV_DATA *V, WENO_DATA *W, PARA *par, int **beta ) {
/*
 * V    - cell averages
 * W    - special representation of cell averages for WENO reconstruction
 * par  - parameter
 * beta - troubled cell indicator
 *
 * save cell averages in an apropriate form for WENO reconstruction to be able to use matrix multiplication for each stencil
 */

    ASSERT ( V );
    ASSERT ( W );
    ASSERT ( par );
    ASSERT ( beta );

    int k, i, s, p;
    int bound;

    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i ) {    // loop cells
            if ( beta[k][i] != OK ) {           // only for troubled cells

                for ( s = 0; s < par->n_stencil; ++s )  // loop stencils
                    for ( p = 0; p < par->Mp1; ++p ) {      // loop cells in stencil

                        // use three or four stencils --> shift difference not equal
                        bound = i - par->M + p + W->shift[s];

                        // treat cells which are not on the 'real' grid
                        if ( strcmp ( par->boundary, "symmetric" ) == 0 || (strcmp ( par->boundary, "hack" ) == 0) ) {
                            if ( bound < 0 ) {
                                gsl_vector_set( W->gslData[k][i][s], p, (double)(par->sym_flags[k]) * V->Avg[k][-1-bound] );
                            } else if (bound > par->N-1)
                                gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][par->N-1] );
                            else
                                gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][bound] );
                        } else if ( strcmp ( par->boundary, "periodic" ) == 0 ) {
                            if ( (bound < 0) || (bound > par->N-1) )
                                bound = translate_boundary ( bound, par->N );

                            ASSERT ( bound >= 0 );
                            ASSERT ( bound < par->N );

                            gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][bound] );
                        } else {
                            if (bound < 0)
                                gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][0] );
                            else if  (bound > par->N-1)
                                gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][par->N-1] );
                            else
                                gsl_vector_set( W->gslData[k][i][s], p, V->Avg[k][bound] );
                        }
                    }
            }
        }
}

// #### #### ####

void Set_Beta_OK ( int **beta, PARA *par ) {
/*
 * beta - troubled cell indicator
 * par  - parameter
 *
 * set all cells to 'OK' except for boundaries
 */

/// ALL EXECPT BOUNDARIES
//    int i,k;
//    for ( k = 0; k < par->var; ++k ) {  // loop variables
//        for ( i = 2; i < par->N-1; ++i )    // loop cells
//            beta[k][i] = OK;

//        beta[k][0] = NEED_AW;
//        beta[k][1] = AW_NEIGH;
//        beta[k][par->N-1] = AW_NEIGH;

//        beta[k][par->N-1] = NEED_AW;
//        beta[k][par->N-2] = AW_NEIGH;
//    }

/// REALLY ALL
    int i, k;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i )      // loop cells
            beta[k][i] = OK;

    OUTPUT_EVO ( "reset beta" );

}

void Set_Beta_Bad ( int **beta, PARA *par ) {
/*
 * beta - troubled cell indicator
 * par  - parameter
 *
 * set all cells troubled
 */

    int i, k;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i )      // loop cells
            beta[k][i] = NEED_AW;

    OUTPUT_EVO ( "set beta bad" );

}

void Set_Beta_Prev ( int **beta, int **beta_prev, PARA *par ) {
/*
 * beta - troubled cell indicator
 * par  - parameter
 *
 * reset previous indicators
 */

    int i, k;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i )      // loop cells
            beta_prev[k][i] = beta[k][i];

    OUTPUT_EVO ( "set previous beta" );

}

void Set_Beta_Sub ( int **beta, int **beta_sub, PARA *par ) {
/*
 * beta - troubled cell indicator on main grid
 * beta_sub - troubled cell indicator of subcells
 * par  - parameter
 *
 * set subcell indicators according to cells on main grid
 */

    int i, k;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->NSC*par->N; ++i ) // loop all subcells
            beta_sub[k][i] = beta[k][(int)(i/par->NSC)];

}

// ### ### ### ### ### ### ###

void Adapt_Troubled_Cells ( int **beta, PARA *par ) {
/*
 * beta - troubled cell indicator
 * par  - parameter
 *
 * if one variable is troubled -> set all other variables for that cell troubled
 */

    int i, k, s;

    for ( i = 1; i < par->N-1; ++i )  // loop cells
        for ( k = 0; k < par->var; ++k )    // loop variables
            if ( beta[k][i] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {  // loop variables
                    beta[s][i] = NEED_AW;

                    if ( beta[s][i-1] != NEED_AW )
                        beta[s][i-1] = AW_NEIGH;
                    if ( beta[s][i+1] != NEED_AW )
                        beta[s][i+1] = AW_NEIGH;
                }

    // boundaries
    if ( strcmp ( par->boundary, "periodic" ) == 0 ) {
        for ( k = 0; k < par->var; ++k ) {
            // i = 0
            if ( beta[k][0] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {
                    beta[s][0] = NEED_AW;

                    if ( beta[s][par->N-1] != NEED_AW )
                        beta[s][par->N-1] = AW_NEIGH;
                }

            // i = N-1
            if ( beta[k][par->N-1] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {
                    beta[s][par->N-1] = NEED_AW;

                    if ( beta[s][0] != NEED_AW )
                        beta[s][0] = AW_NEIGH;
                }
        }
    }
    else if ( (strcmp ( par->boundary, "symmetric" ) == 0) || (strcmp ( par->boundary, "none" ) == 0 ) ) {
        for ( k = 0; k < par->var; ++k ) {
            // i = 0
            if ( beta[k][0] == NEED_AW )
                for ( s = 0; s < par->var; ++s )
                    beta[s][0] = NEED_AW;


            // i = N-1
            if ( beta[k][par->N-1] == NEED_AW )
                for ( s = 0; s < par->var; ++s )
                    beta[s][par->N-1] = NEED_AW;
        }
    }
    else {
        printf ( "*** ERROR - specify boundary condition!\n" );
        ASSERT ( 0 );
    }

}

void Adapt_Troubled_Cells_extended ( int **beta, PARA *par ) {
/*
 * beta - troubled cell indicator
 * par  - parameter
 *
 * if one variable is troubled -> set all other variables for that cell troubled
 * set neighbours also troubled and create new neighbours...
 */

    int i, k, s;

    for ( i = 1; i < par->N-1; ++i )  // loop cells
        for ( k = 0; k < par->var; ++k )    // loop variables
            if ( beta[k][i] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {  // loop variables
                    beta[s][i] = NEED_AW;

                    if ( beta[s][i-1] != NEED_AW )
                        beta[s][i-1] = AW_NEIGH;

                    if ( beta[s][i+1] != NEED_AW )
                        beta[s][i+1] = AW_NEIGH;
                }

    // set 'old' neighbours troubled
    for ( i = 0; i < par->N; ++i )
        if ( beta[0][i] == AW_NEIGH )
            beta[0][i] = NEED_AW;

    // find new neighbours
    for ( i = 1; i < par->N-1; ++i )  // loop cells
        if ( beta[0][i] == NEED_AW )
            for ( s = 0; s < par->var; ++s ) {  // loop variables
                beta[s][i] = NEED_AW;

                if ( beta[s][i-1] != NEED_AW )
                    beta[s][i-1] = AW_NEIGH;

                if ( beta[s][i+1] != NEED_AW )
                    beta[s][i+1] = AW_NEIGH;
            }

    // boundaries
    if ( strcmp ( par->boundary, "periodic" ) == 0 ) {
        for ( k = 0; k < par->var; ++k ) {
            // i = 0
            if ( beta[k][0] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {
                    beta[s][0] = NEED_AW;

                    if ( beta[s][par->N-1] != NEED_AW )
                        beta[s][par->N-1] = AW_NEIGH;
                }

            // i = N-1
            if ( beta[k][par->N-1] == NEED_AW )
                for ( s = 0; s < par->var; ++s ) {
                    beta[s][par->N-1] = NEED_AW;

                    if ( beta[s][0] != NEED_AW )
                        beta[s][0] = AW_NEIGH;
                }
        }
    }
    else if ( (strcmp ( par->boundary, "symmetric" ) == 0) || (strcmp ( par->boundary, "none" ) == 0 ) ) {
        for ( k = 0; k < par->var; ++k ) {
            // i = 0
            if ( beta[k][0] == NEED_AW )
                for ( s = 0; s < par->var; ++s )
                    beta[s][0] = NEED_AW;


            // i = N-1
            if ( beta[k][par->N-1] == NEED_AW )
                for ( s = 0; s < par->var; ++s )
                    beta[s][par->N-1] = NEED_AW;
        }
    }
    else {
        printf ( "*** ERROR - specify boundary condition!\n" );
        ASSERT ( 0 );
    }

}

void Count_Troubled_Cells ( int *beta, PARA *par ) {
/*
 * beta - troubled cell indicator of first variable
 * par - parameter
 *
 * count number of troubled cells and add to average in par
 */

    int i;
    int temp_t = 0;
    int temp_n = 0;

    for ( i = 0; i < par->N; ++i ) {
        if ( beta[i] == NEED_AW )
            temp_t++;
        else if ( beta[i] == AW_NEIGH )
            temp_n++;
    }

    par->N_trouble += (double)temp_t/((double)par->N);
    par->N_neigh   += (double)temp_n/((double)par->N);
}







