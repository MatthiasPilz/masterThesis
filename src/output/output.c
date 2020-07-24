#include "ader.h"
#include "output.h"



// ### ### ### ### ### ### ###

void Output_Parameter ( PARA *par ) {

    printf ( "\n\tSimulation of %s Equation using the %s-DG scheme.\n", par->project, par->scheme );
	
    printf ( "\n\t#######--- Parameter ---##########\n" );
    printf ( "\t#\tM\t= %8d\n", par->M);
    printf ( "\t#\tproject\t= %s\n", par->project );
    printf ( "\t#\tshape\t= %s\n", par->shape );
    printf ( "\t#\n" );

    printf ( "\t# GRID:\n" );
    printf ( "\t#\tN\t= %8d\n", par->N);
    printf ( "\t#\txRange\t= %6lf\n", par->xRange );
    printf ( "\t#\txMin\t= %6lf\n", par->xMin );
    printf ( "\t#\txMax\t= %6lf\n", par->xMax );
    printf ( "\t#\tdx\t= %6lf\n", par->dx );
    printf ( "\t#\tCFL\t= %6lf\n", par->CFL );
    printf ( "\t#\tboundary\t= %s\n", par->boundary );
    printf ( "\t#\n" );

    printf ( "\t# TIME:\n" );
    printf ( "\t#\tNt\t= %8d\n", par->Nt);
    printf ( "\t#\ttEND\t= %6lf\n", par->tEND );
    printf ( "\t#\tdt\t= %6lf\n", par->dt );
    printf ( "\t#\n" );

    if ( strcmp ( par->scheme, "dumbser" ) == 0 ) {
        printf ( "\t# ITERATION\n" );
        printf ( "\t#\titerMAX\t= %8d\n", par->iterMAX);
        printf ( "\t#\teps\t= %.5e\n", par->eps);
        printf ( "\t#\n" );
    }
    else if ( strcmp ( par->scheme, "RK" ) == 0 )
        printf ( "\t#\tstages\t= %d\n", par->stages );

    printf ( "\t# TROUBLED CELL INDICATOR:\n" );
    printf ( "\t#\tkappa\t= %6lf\n", par->kappa );
    printf ( "\t#\n" );

#ifdef CONVERGENCE
        printf ( "\t# CONVERGENCE:\n" );
        printf ( "\t#\tM_start\t= %d\n", par->conv_M_start );
        printf ( "\t#\tM_end\t= %d\n", par->conv_M_end );
        printf ( "\t#\tM_step\t= %d\n", par->change_M );
        printf ( "\t#\tN_start\t= %d\n", par->conv_N_start );
        printf ( "\t#\tN_start\t= %d\n", par->conv_N_end );
        printf ( "\t#\tN_factor\t= %d\n", par->change_N );
        printf ( "\t#\n" );
#endif

    printf ( "\t##################################\n\n" );
	
}

void Setup_Output_Data_Structure ( PARA *par ) {

    char command[255];

    // setup data structure
    printf ( "\nSetting up data structure for output:\n" );

    sprintf ( command, "mkdir ../data/%s", par->project );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s\n", par->project );

    sprintf ( command, "mkdir ../data/%s/%s", par->project, par->shape );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s/%s\n", par->project, par->shape );

    sprintf ( command, "mkdir ../data/%s/%s/%s", par->project, par->shape, par->scheme );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s/%s/%s\n", par->project, par->shape, par->scheme );

    sprintf ( command, "../data/%s/%s/%s", par->project, par->shape, par->scheme );
    Copy_Parameter_Data ( par->file, command );

    printf ( "\n" );

}

// ### ### ### ### ### ### ###

static double lagrange_polynom ( int k, double x, const int M, double *LAM ) {

    int i;
    double result = 1.0;
    int order = M+1;

    // product loop
    for ( i = 0; i < order; ++i ) {
        if ( i != k )
            result *= ( (x - LAM[i]) / (LAM[k] - LAM[i]) );
    }

    return result;
}

void Output_Polynom_Solution ( gsl_vector ***W, double *LAM, PARA *par, char *name ) {
    int i, j, k, s;
    double dx_plot = par->dx / (par->Npoly+1.0);   ASSERT ( dx_plot > 0.0 );
    double y;
    char filename[255];

    for ( s = 0; s < par->var; ++s ) {
        sprintf ( filename, "../data/%s/%s/%s/%d_poly%d_%s.txt", par->project, par->shape, par->scheme, par->N, s, name );

        FILE *fp = fopen ( filename, "a" );
        ASSERT ( fp );

        fprintf ( fp, "\"Time=%lf\n", par->time );

        for ( i = 0; i < par->N; ++i )
            for ( j = 1; j < par->Npoly+1; ++j ) {
                y = 0.0;
                ASSERT ( j/(par->Npoly+1.0) >= 0.0 );
                ASSERT ( j/(par->Npoly+1.0) <= 1.0 );

                for ( k = 0; k < par->Mp1; ++k )
                    y += ( lagrange_polynom ( k, j/(par->Npoly+1.0), par->M, LAM ) * gsl_vector_get ( W[s][i], k ) );

                fprintf ( fp, "%.15f %.15f\n", par->xMin + i*par->dx+j*dx_plot, y );
            }

        fclose ( fp );
    }

// 	printf ( "# SUCCESS # - wrote polynomial reconstruction in file: %s\n", filename );

}

void Output_Polynom_Solution_on_Lagrange_Points ( gsl_vector ***W, double *LAM, PARA *par, char *name ) {

    int i, j, k;
    char filename[255];

    for ( k = 0; k < par->var; ++k ) {
        sprintf ( filename, "../data/%s/%s/%s/%d_LAM%d_%s.txt", par->project, par->shape, par->scheme, par->N, k, name );

        FILE *fp = fopen ( filename, "a" );
        ASSERT ( fp );

        fprintf ( fp, "\"Time=%lf\n", par->time );

        for ( i = 0; i < par->N; ++i ) {

            for ( j = 0; j < par->Mp1; ++j ) {

                fprintf ( fp, "%.15f %15f\n", par->xMin + (i+LAM[j])*par->dx, gsl_vector_get ( W[k][i], j ) );

            }

        }

    }

}

// ### ### ### ### ### ### ###

void Output_GSL_Matrix ( gsl_matrix *A, int sizeR, int sizeC ) {

    ASSERT ( A );
    ASSERT ( sizeR > 0 );
    ASSERT ( sizeC > 0 );

    int i, j;
    for ( i = 0; i < sizeR; ++i ) {
        for ( j = 0; j < sizeC; ++j )
            printf ( "%.4f\t", gsl_matrix_get ( A, i, j ) );

        NEWLINE
    }

}

// ### ### ### ### ### ### ###

void Output_Double_Matrix ( double **A, int sizeR, int sizeC ) {

    ASSERT ( A );
    ASSERT ( sizeR > 0 );
    ASSERT ( sizeC > 0 );

    int i, j;
    for ( i = 0; i < sizeR; ++i ) {
        for ( j = 0; j < sizeC; ++j )
            printf ( "%.20f\t", A[i][j] );

        NEWLINE
    }

}

// ### ### ### ### ### ### ###

void Output_Evolution_Summary ( PARA *par ) {

    ASSERT ( par );
    NEWLINE

    par->chanAvg /= (double)par->Nt;
    par->iterAvg /= (double)par->Nt;

    if ( par->iterFail )
        printf ( "Iteration FAILED - for at least one cell\n" );

    printf ( "average number of iterations = %.5e\n", par->iterAvg );
    printf ( "average change accuracy = %.5e\n", par->chanAvg );
    printf ( "maximum change = %.5e\n", par->changeMAX );

}

// ### ### ### ### ### ### ###

void Output_Current_Status ( PARA *par ) {

    ASSERT ( par );
    NEWLINE

    printf ( "current time: t = %.5lf\n", par->time );
    printf ( "%d%% of simulation done!\n", (int)( 100*par->count/par->Nt ) );

    if ( strcmp ( par->scheme, "dumbser" ) == 0 ) {
//        printf ( "maximum change so far = %.5e\n", par->changeMAX );
        printf ( "average change so far = %.5e\n", par->chanAvg/par->count );
        printf ( "average number of iterations = %.5f\n", par->iterAvg/par->count );
    }

    printf ( "average number of troubled cells = %.5f\n", par->N_trouble/((double)par->count) );
    printf ( "average number of neighbour cells = %.5f\n", par->N_neigh/((double)par->count) );

}

// ### ### ### ### ### ### ###

void Output_Current_Status_Subcell ( PARA *par ) {

    ASSERT ( par );
    NEWLINE

    printf ( "current time: t = %.5lf\n", par->time );
    printf ( "%d%% of simulation done!\n", (int)( 100*par->count/par->Nt ) );

    if ( strcmp ( par->scheme, "dumbser" ) == 0 ) {
//        printf ( "maximum change so far = %.5e\n", par->changeMAX );
        printf ( "average change so far = %.5e\n", par->chanAvg/par->count );
        printf ( "average number of iterations = %.5f\n", par->iterAvg/(par->N_trouble + par->N_neigh) );
    }

    printf ( "average number of troubled cells = %.5f\n", par->N_trouble/((double)par->count) );
    printf ( "average number of neighbour cells = %.5f\n", par->N_neigh/((double)par->count) );

}

// ### ### ### ### ### ### ###

void Output_Final_Status ( PARA *par, PARA *par_AW ) {

    ASSERT ( par );
    NEWLINE
    NEWLINE
    printf ( "FINAL TIME reached:\n" );

    NEWLINE NEWLINE
    printf ( "GRID status:\n" );
    Output_Current_Status ( par );
    NEWLINE
    printf ( "total number of evals: %ld\n", par->evalCountF + par->evalCountS );
    printf ( "total number of eval per cell: %.2f\n", (double)(par->evalCountF+par->evalCountS)/par->N );
    NEWLINE
    printf ( "(predictor) number of func_evals: %ld\n", par->evalPrediF );
    printf ( "(predictor) number of f_eval per cell: %.2f\n", (double)par->evalPrediF/par->N );
    NEWLINE
    printf ( "number of func_evals per cell per step: %.2f\n", (double)par->evalCountF/(par->N*par->count) );
    printf ( "(predictor) number of f_eval per cell per step: %.2f\n", (double)par->evalPrediF/(par->N*par->count) );

    if ( par->source_flag ) {
        NEWLINE
        printf ( "(predictor) number of source_evals: %ld\n", par->evalPrediS );
        printf ( "(predictor) number of s_eval per cell: %.2f\n", (double)par->evalPrediS/par->N );
        NEWLINE
        printf ( "number of source_evals per cell per step: %.2f\n", (double)par->evalCountS/(par->N*par->count) );
        printf ( "(predictor) number of s_eval per cell per step: %.2f\n", (double)par->evalPrediS/(par->N*par->count) );
    }



    NEWLINE
    NEWLINE



    printf ( "SUBCELL status:\n" );
    Output_Current_Status_Subcell ( par_AW );
    NEWLINE
    printf ( "total number of evals: %ld\n", par_AW->evalCountF + par_AW->evalCountS );
    printf ( "total number of eval per troubled cell: %.2f\n", (double)((par_AW->evalCountF+par->evalCountS)*par_AW->count)/(par_AW->N*(par->N_trouble + par->N_neigh)) );
    NEWLINE
    printf ( "(predictor) number of func_evals: %ld\n", par_AW->evalPrediF );
    printf ( "(predictor) number of f_eval per troubled cell: %.2f\n", (double)(par_AW->evalPrediF*par_AW->count)/(par_AW->N*(par->N_trouble + par->N_neigh)) );
    NEWLINE
    printf ( "number of func_evals per troubled cell per step: %.2f\n", (double)par_AW->evalCountF/(par_AW->N*(par->N_trouble + par->N_neigh)) );
    printf ( "(predictor) number of f_eval per troubled cell per step: %.2f\n", (double)par_AW->evalPrediF/(par_AW->N*(par->N_trouble + par->N_neigh)) );

    if ( par_AW->source_flag ) {
        NEWLINE
        printf ( "(predictor) number of source_evals: %ld\n", par_AW->evalPrediS );
        printf ( "(predictor) number of s_eval per troubled cell: %.2f\n", (double)par_AW->evalPrediS/(par_AW->N*(par->N_trouble + par->N_neigh)) );
        NEWLINE
        printf ( "number of source_evals per troubled cell per step: %.2f\n", (double)par_AW->evalCountS/(par_AW->N*(par->N_trouble + par->N_neigh)) );
        printf ( "(predictor) number of s_eval per troubled cell per step: %.2f\n", (double)par_AW->evalPrediS/(par_AW->N*(par->N_trouble + par->N_neigh)) );
    }
}

// ### ### ### ### ### ### ###

void Output_Troubled_Cells ( int **beta, PARA *par, char *name ) {

    int i, j;
    char filename[255];

    sprintf ( filename, "../data/%s/%s/%s/%d_beta%d_%s.txt", par->project, par->shape, par->scheme, par->N, j, name );

    FILE *fp = fopen ( filename, "a" );

    fprintf ( fp, "\"Time=%lf\n", par->time );

    for ( i = 0; i < par->N; ++i )
       fprintf ( fp, "%.15f %.3f\n", par->xMin + (i+0.5)*par->dx, 0.5*beta[0][i] );

    fclose ( fp );
}

// ### ### ### ### ### ### ###

void Output_Troubled_Cells_Conv ( int **beta, PARA *par, char *name ) {

    int i;
    char filename[511];

    sprintf ( filename, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/%s.txt", par->project, par->M_sub, par->shape, par->scheme, (2*par->M+1)*par->CFL, (2*par->M+1)*par->N, name );

    FILE *fp = fopen ( filename, "a" );

    fprintf ( fp, "\"Time=%lf\n", par->time );

    for ( i = 0; i < par->N; ++i )
       fprintf ( fp, "%.15f %.3f\n", par->xMin+(i+0.5)*par->dx, 0.5*beta[0][i] );

    fclose ( fp );
}

// ### ### ### ### ### ### ###

void Output_Troubled_Cells_QS ( int **beta, PARA *par, char *name ) {
/*
 * Output for troubled cell plots as in Qiu and Shu 2005
 * without time stamp and not at t=0
 */

    int i, j;
    char filename[255];

    if ( par->time != 0.0 ) {

    for ( j = 0; j < par->var; ++j ) {
        sprintf ( filename, "../data/%s/%s/%s/%d_beta%d_%s.txt", par->project, par->shape, par->scheme, par->N, j, name );

        FILE *fp = fopen ( filename, "a" );

        for ( i = 0; i < par->N; ++i ) {
            if ( beta[j][i] == NEED_AW ) {
                fprintf ( fp, "%.15f %.3f\n", par->xMin + (i+0.5)*par->dx, par->time*0.5*beta[j][i] );
            }
        }
        fclose ( fp );
    }
    }


}












