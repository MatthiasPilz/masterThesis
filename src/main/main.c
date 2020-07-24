/*

  ADER-DG algorithm with a posteriori subcell limiter to solve hyperbolic balance law equations

  different equations are included as projects

  complete theoretical background is given in papers by Dumbser et al.
  [1] - 'ADER-WENO Finite Volume Scheme with Space-Time Adaptive Mesh Refinement',2013
  [2] - 'A Posteriori Subcell Limiting of the DG FEM for Hyperbolic Conservation Laws',2014

*/


#include "ader.h"
#include "main.h"


int main ( int argc, char *argv[] ) {

/// DEFINITIONS
// set DG parameter
    PARA par;
    Read_Parameter                  ( argv[1], &par );

// set function pointer
    PROJ pro = Set_Function_Pointer ( &par );
    pro.f_var_num                   ( &par );
    Set_Scheme_Function_Pointer     ( &par );

// set ADER WENO parameter
    PARA par_AW;
    Delete_old_Iteration_File       ( &par );


#ifdef CONVERGENCE
    char output_name[255];

  for ( par.M = par.conv_M_start; par.M <= par.conv_M_end; par.M += par.change_M ) {
    for ( par.N = par.conv_N_start; par.N <= par.conv_N_end; par.N *= par.change_N ) {
        Set_Convergence_Parameter ( &par );

#endif

    if ( strcmp( par.boundary, "hack" ) == 0 )
        Hack_Set_DG_Parameter   ( &par );

    Set_ADER_WENO_Parameter ( &par, &par_AW );
//    Set_RK_Dense_Output_Function ( &par );
    Output_Parameter ( &par_AW );

#ifndef CONVERGENCE
            // prepare file structure
            Delete_Data ( &par );
            Setup_Output_Data_Structure ( &par );
#endif


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


        /// free regular grid
        /// adjust parameter struct (N, xMax, xRange, .?. )
        /// continue regular treatment


    if ( strcmp( par.boundary, "hack" ) == 0 )
        Hack_Save_TOV_Boundary_Flux ( V, W, A, DG, &par, &pro );



#ifdef CONVERGENCE
    Setup_Convergence_Data_Structure    ( &par_AW );
    pro.f_convergence_output            ( V_AW, &par_AW, "ana" );
    printf ( "NT = %d\n", par.Nt );
#endif


/// TIME EVOLUTION
    for ( par.count = 0; par.count < par.Nt; par.count++ ) {
            par.evalTempF = 0;
            par.evalTempS = 0;
        Update_Time_Parameter ( &par, &par_AW );
        // output status
        if ( !(par.count % ((int)(10*par_AW.N)) ) ) {
            Output_Current_Status ( &par );
        }


        // interpolation of current data to subcells and cell average on grid
        Interpolation               ( DG->U, V_AW->Avg_adm, DG->Int, &par, &par_AW );
        Interpolation_Good_Cells    ( DG, V_AW->Avg, DG->Int, &par, &par_AW );
        Cell_Averages               ( DG->U, DG->V_grid, DG->AVG, &par );


    /// OUTPUT
        if ( !( par.count % par.output ) ) {
#ifndef CONVERGENCE
            pro.f_output_averages   ( V_AW->Avg, &par_AW, "avg" );
            pro.f_output_averages   ( DG->V_grid, &par, "U_avg" );
//            Output_Polynom_Solution ( DG->U, DG->LAM, &par, "U" );
//            Output_Troubled_Cells   ( DG->beta_sub, &par_AW, "beta_sub" );
            Output_Troubled_Cells   ( DG->beta, &par, "beta" );

#else
            sprintf ( output_name, "Nt_%08d", par.count );
            pro.f_convergence_output    ( V_AW, &par_AW, output_name );
            Output_Troubled_Cells_Conv  ( DG->beta, &par, "beta" );
#endif
        }


        Set_Beta_Bad        ( DG->beta, &par );
        // local DG-predictor (either Iteration or RK)
        f_predictor_step    ( A, DG->U, DG->beta, &par, &pro );
            par.evalPrediF += par.evalTempF;
            par.evalPrediS += par.evalTempS;

        // ADER step
        ADER_DG_Step        ( DG, A, &par, &pro );

        // interpolation of predictor data to subcells and averages
        Interpolation       ( DG->pU, DG->V, DG->Int, &par, &par_AW );
        Cell_Averages       ( DG->U,  DG->V_grid_pre, DG->AVG, &par );

        // admissibility test (DUMBSER  or  COCKBURN & SHU)
        Set_Beta_OK                     ( DG->beta, &par );
//        Numerical_Admissibility         ( DG, V_AW->Avg_adm, &par_AW );   /// DUMBSER
        Troubled_Cell_Indicator         ( DG, &par );                     /// COCKBURN

        if ( pro.f_phy_admiss != NULL )
            pro.f_phy_admiss ( DG, &par );

            //Output_Troubled_Cells_QS  ( DG->beta, &par, "troubled_cells" );

//        Set_Beta_Bad                    ( DG->beta, &par );   /// uncomment for setting all cells troubled
//        Set_Beta_OK                     ( DG->beta, &par );   /// uncomment for setting no cells troubled

        Adapt_Troubled_Cells_extended   ( DG->beta, &par );
        Count_Troubled_Cells            ( DG->beta[0], &par );

        // accept solution?
        Accept_Unlimited_DG ( DG, V_AW, W_AW, A_AW, &par, &par_AW, &pro );
        Set_Beta_Prev       ( DG->beta, DG->beta_prev, &par );

            par.evalCountF += par.evalTempF;
            par.evalCountS += par.evalTempS;

    }


/// FINAL TIME STEP TO REACH EXACT tEND
    par.time = par.dt * par.count;
    par_AW.time = par.time;
    Interpolation_Good_Cells    ( DG, V_AW->Avg, DG->Int, &par, &par_AW );

    if ( par.time < par.tEND ) {
        Final_Step                  ( V_AW, W_AW, A_AW, DG->beta_sub, &par_AW, &pro );
    }
    par_AW.N_trouble = par.N_trouble;
    par_AW.N_neigh   = par.N_neigh;
    Output_Final_Status ( &par, &par_AW );

#ifdef CONVERGENCE
        Save_Iteration_Statistic ( &par, &par_AW );
        pro.f_convergence_output ( V_AW, &par_AW, "final" );
        NEWLINE NEWLINE NEWLINE NEWLINE NEWLINE
#else
        par.time = par.tEND;
        pro.f_output_averages   ( V_AW->Avg, &par_AW, "avg_final" );
#endif
    timer_print ( &par );
//    Save_Iteration_Statistic ( &par, &par_AW );

/// FREE DATA
//    NEWLINE
//    Free_FV_DATA ( V, &par );
//    Free_WENO_DATA ( W, &par );
//    Free_ADER_DATA ( A, &par );
//    Free_DG_DATA ( DG, &par );

//    Free_FV_DATA ( V_AW, &par_AW );
//    Free_WENO_DATA ( W_AW, &par_AW );
//    Free_ADER_DATA ( A_AW, &par_AW );



#ifdef CONVERGENCE
    }
  }

  Calc_Convergence_Test ( par.file );
#endif
    return 0;
}





/// WRAPPER

int SUBCELL_WENO_TIME_STEP ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro ) {

    par->evalTempF = 0;
    par->evalTempS = 0;

    // WENO reconstruction
    FV_2_WENO_DATA          ( V, W, par, beta );
    Solve_Algebraic_System  ( W, par, beta );
    Calculate_WENO_Coeff    ( W, par, beta );

    // local predictor
    f_predictor_step        ( A, W->WT, beta, par, pro );
    par->evalPrediF += par->evalTempF;
    par->evalPrediS += par->evalTempS;

    // FV scheme to update data
    Calculate_Flux_Averages ( V, A, par, beta, pro );

    if ( par->source_flag )
        Calculate_Source_Averages ( V, A, par, beta, pro );

    Update_Finite_Volume_Scheme ( V, par, beta );
    par->evalCountF += par->evalTempF;
    par->evalCountS += par->evalTempS;

    return par->iterFail;
}


int ALLOC_INIT_ALL_DATA ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, DG_DATA *DG, PARA *par, PROJ *pro ) {

// create initial data, allocate and initialize FV_DATA
    Allocate_FV_DATA_general ( V, par );
    Init_FV_DATA_general ( V, par );
    pro->f_init_data ( V, par );
    NEWLINE

// allocate and initialize WENO_DATA
    Allocate_WENO_DATA_general ( W, par );
    Init_WENO_DATA_general ( W, par );
    NEWLINE

// allocate and initialize ADER_DATA
    Allocate_ADER_DATA_general ( A, par );
    Init_ADER_DATA_general ( A, par );
    NEWLINE

// allocate and initialize DG_DATA
    Allocate_DG_DATA_general ( DG, par );
    Init_DG_DATA_general ( V, W, DG, par, pro );

    return 0;

}


int ALLOC_INIT_ADER_WENO ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, PARA *par, PROJ *pro ) {

// create initial data, allocate and initialize FV_DATA
    Allocate_FV_DATA_general ( V, par );
    Init_FV_DATA_general ( V, par );
    pro->f_init_data ( V, par );
    NEWLINE

// allocate and initialize WENO_DATA
    Allocate_WENO_DATA_general ( W, par );
    Init_WENO_DATA_general ( W, par );
    NEWLINE

// allocate and initialize ADER_DATA
    Allocate_ADER_DATA_general ( A, par );
    Init_ADER_DATA_general ( A, par );

    return 0;

}

int SUBCELL_WENO_TIME_STEP_woStatistics ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro ) {

    // WENO reconstruction
    FV_2_WENO_DATA          ( V, W, par, beta );
    Solve_Algebraic_System  ( W, par, beta );
    Calculate_WENO_Coeff    ( W, par, beta );

    // local predictor
    f_predictor_step        ( A, W->WT, beta, par, pro );

    // FV scheme to update data
    Calculate_Flux_Averages ( V, A, par, beta, pro );

    if ( par->source_flag )
        Calculate_Source_Averages ( V, A, par, beta, pro );

    Update_Finite_Volume_Scheme ( V, par, beta );

    return par->iterFail;
}















