// code snippets needed for testing --> do not compile!


/// #### #### #### CONVERGENCE TEST #### #### ####

// definitions
#ifdef CONVERGENCE
  //  for ( par.CFL = 1.0; par.CFL >= 0.01; par.CFL /= 2 ) {
        for ( par.N = par.conv_N_start; par.N <= par.conv_N_end; par.N += par.conv_step ) {
            if ( par.N <= 1000 && par.N >= 90 )
                par.conv_step = 90;
            if ( par.N > 1000 )
                par.conv_step = 900;

            par.dx = (double) par.xRange / par.N;
            par.dt = ( par.CFL * 1.0 / par.N );
            par.Nt = par.tEND / par.dt;

            par.chanAvg = 0.0;
            par.iterAvg = 0.0;
            par.iterFail = 0;
            par.changeMAX = 0.0;
#endif
// allocation
#ifdef CONVERGENCE
    Setup_Convergence_Data_Structure ( &par );
    pro.f_convergence_output ( V, &par, "init" );
    printf ( "NT = %d\n", par.Nt );
#endif
// evolution
#ifdef CONVERGENCE
    pro.f_convergence_output ( V, &par, "final" );
#endif
// free
#ifdef CONVERGENCE
        }
  //  }
#endif

/// #### #### #### #### #### #### #### #### #### ####




/// TESTS

double gauss ( double x ) {
    return ( 2.0*exp ( -pow((x - 0.5),2) * 150.0 ) );
}

if ( strcmp ( par.project, "burger" ) == 0 ) {
    /// TEST - Burger!
    int i;
    FILE *fp1 = fopen ( "../data/Test.txt", "w" );
    FILE *fp2 = fopen ( "../data/Test_vgl.txt", "w" );
    for ( i = 0; i < par.N; ++i )
        fprintf ( fp1, "%lf %lf\n", par.xMin + (i+0.5)*par.dx - V->Avg[0][i]*par.tEND , gauss ( par.xMin + (i+0.5)*par.dx - V->Avg[0][i]*par.tEND ) );

    for ( i = 0; i < par.N * 10; ++i )
        fprintf ( fp2, "%lf %lf\n", par.xMin + (i+0.5)*par.dx/10.0, gauss (par.xMin +  (i+0.5)*par.dx/10.0) );

    fclose ( fp1 );
    fclose ( fp2 );
}
else if ( strcmp ( par.project, "dumbser" ) == 0 ) {

    /// TEST - Dumbser!
    FILE *fp1 = fopen ( "../data/Test.txt", "w" );
    int i;
    for ( i = 0; i < par.N*10; ++i )
        fprintf ( fp1, "%lf %lf\n", (i+0.5)*par.dx/10.0, 0.1*Dumbser_Func((i+0.5)*par.dx/10.0, 0.119, 0));

    fclose ( fp1 );

}




void TEST_WENO_gslData_Output ( WENO_DATA *W, PARA *par ) {

    int i, s, p;
    char name[255];
    int bound;

    for ( i = 0; i < par->N; ++i )
        for ( s = 0; s < par->n_stencil; ++s ) {

            sprintf ( name, "../data/test/M%d_i%d_s%d.txt", par->M, i, s );
            FILE *fp = fopen ( name, "w" );
            for ( p = 0; p < par->Mp1; ++p ) {
                bound = i - par->M + s + p;
                if ( (bound < 0) || (bound > par->N-1) )
                        bound = translate_boundary ( bound, par->M, par->N );

                ASSERT ( bound >= 0 );
                ASSERT ( bound < par->N );

                fprintf ( fp, "%.15f %.15f\n", (bound+0.5)*par->dx, gsl_vector_get ( W->gslData[0][i][s], p  ) );
            }
            fclose ( fp );

        }

}











/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// komplettes ADER_WENO Programm mit Konvergenztest

int main ( int argc, char *argv[] ) {

/// DEFINITIONS
// set DG parameter
    PARA par;
    Read_Parameter ( argv[1], &par );

// set function pointer
    PROJ pro = Set_Function_Pointer ( &par );
    pro.f_var_num ( &par );

// prepare file structure
    Delete_Data ();
    Setup_Output_Data_Structure ( &par );


#ifdef CONVERGENCE
  //  for ( par.CFL = 1.0; par.CFL >= 0.01; par.CFL /= 2 ) {
        for ( par.N = par.conv_N_start; par.N <= par.conv_N_end; par.N += par.conv_step ) {
            if ( par.N <= 112 ) { par.conv_step = 16; }
            else if ( par.N > 112 && par.N < 512 ) { par.conv_step = 64; }
            else if ( par.N >= 512 ) { par.conv_step = 128; }
            else { printf ( "ERROR \n" ); }

            par.dx = (double) par.xRange / par.N;
            par.dt = ( par.CFL * 1.0 / par.N );
            par.Nt = par.tEND / par.dt;

            par.chanAvg = 0.0;
            par.iterAvg = 0.0;
            par.iterFail = 0;
            par.changeMAX = 0.0;
#endif

/// ALLOCATION AND INITIALIZATION
    FV_DATA     *V  = malloc ( sizeof ( FV_DATA ) );
    WENO_DATA   *W  = malloc ( sizeof ( WENO_DATA ) );
    ADER_DATA   *A  = malloc ( sizeof ( ADER_DATA ) );
    ALLOC_INIT_ADER_WENO ( V, W, A, &par, &pro );

    NEWLINE

#ifdef CONVERGENCE
    Setup_Convergence_Data_Structure ( &par );
    pro.f_convergence_output ( V, &par, "init" );
    printf ( "NT = %d\n", par.Nt );
#endif

/// TIME EVOLUTION
    ADER_WENO_EVOLUTION ( V, W, A, &par, &pro );

#ifdef CONVERGENCE
    pro.f_convergence_output ( V, &par, "final" );
#endif

/// FREE DATA
    free ( V );
    free ( W );
    free ( A );

#ifdef CONVERGENCE
        }
  //  }
#endif

    return 0;

}



















// ADER DG ANFANG
int main ( int argc, char *argv[] ) {

/// DEFINITIONS
// set DG parameter
    PARA par;
    Read_Parameter ( argv[1], &par );

// set function pointer
    PROJ pro = Set_Function_Pointer ( &par );
    pro.f_var_num ( &par );

// set ADER WENO parameter
    PARA par_AW;
    Set_ADER_WENO_Parameter ( &par, &par_AW );
    Output_Parameter ( &par_AW );


// prepare file structure
    Delete_Data ();
    Setup_Output_Data_Structure ( &par );


/// ALLOCATION AND INITIALIZATION
    FV_DATA     *V  = malloc ( sizeof ( FV_DATA ) );
    WENO_DATA   *W  = malloc ( sizeof ( WENO_DATA ) );
    ADER_DATA   *A  = malloc ( sizeof ( ADER_DATA ) );
    DG_DATA     *DG = malloc ( sizeof ( DG_DATA ) );
    ALLOC_INIT_ALL_DATA ( V, W, A, DG, &par, &pro );

    NEWLINE

    FV_DATA     *V_AW  = malloc ( sizeof ( FV_DATA ) );
    WENO_DATA   *W_AW  = malloc ( sizeof ( WENO_DATA ) );
    ADER_DATA   *A_AW  = malloc ( sizeof ( ADER_DATA ) );
    ALLOC_INIT_ADER_WENO ( V_AW, W_AW, A_AW, &par_AW, &pro );


/// TIME EVOLUTION
//    for ( par.count = 0; par.count < par.Nt; par.count++ ) {
        par.time = par.dt * par.count;

        // interpolation of current data to subcells
        Interpolation ( DG->U, V_AW->Avg, &par, &par_AW );

//        Output_Cell_Averages ( V_AW, &par_AW );
        Output_Polynom_Solution ( DG->U, &par, 1 );

        // local DG-predictor + ADER step
        pro.f_iteration     ( A, DG->U, &par );
        pro.f_ader_dg_step  ( DG, A, &par );

        par.time += 0.1;
        Output_Polynom_Solution ( DG->U, &par, 1 );
        par.time += 0.1;
        Output_Polynom_Solution ( DG->pU, &par, 1 );

        // admissibility test

//    }


//    ADER_WENO_EVOLUTION ( V_DG, W_DG, A_DG, &par_DG, &pro );



/// FREE DATA
    free ( V );
    free ( W );
    free ( A );
    free ( DG );


    return 0;

}





void Set_Initial_Iteration_Condition ( ADER_DATA *A, gsl_vector ***W, PARA *par ) {
/*
 * A   - Vold are current coefficients for space-time polynomial
 * W   - coefficients for spatial polynomials from WENO reconstruction
 * par - parameter
 *
 * set space-time coefficients as if we have a steady solution --> all time-'layers' have same data
 */

    ASSERT ( A );
    ASSERT ( W );
    ASSERT ( par );

    int k, i, j;
    for ( k = 0; k < par->var; ++k )    // loop variables
        for ( i = 0; i < par->N; ++i )      // loop cells
            for ( j = 0; j < par->size; ++j )   // loop gaussian quadrature points
                gsl_vector_set ( A->Vold[k][i], j, gsl_vector_get( W[k][i], j%(par->Mp1) ) );

    OUTPUT_EVO ( "set initial iteration conditions" );

}

// ####











        par.time += 0.1;
        Output_Polynom_Solution ( DG->U, &par, 1 );
        par.time += 0.1;
        Output_Polynom_Solution ( DG->pU, &par, 1 );

        int i;
        FILE *fp = fopen ( "../data/beta.txt", "w" );
        fprintf ( fp, "\"Time=%lf\n", par.time );
        for ( i = 0; i < par.N; ++i )
            fprintf ( fp, "%lf %d\n", (i+0.5)*par.dx, DG->beta[0][i] );
        fclose ( fp );


