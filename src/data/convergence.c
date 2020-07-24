#include "ader.h"
#include "data.h"



void Set_Convergence_Parameter ( PARA *par ) {

    par->size = pow ( par->M+1, 2 );
    par->Mp1 = par->M+1;
    if ( par->M % 2 ) par->n_stencil = 4;
    else par->n_stencil = 3;

    if ( par->M == 1 )
        par->n_stencil = 2;

    par->NSC = 2*par->M + 1;
    par->CFL = par->CFL_sub / ( 2.0*par->M + 1 );
    par->iterMAX = par->Mp1 + par->iterEXTRA;
    par->dx = (double) par->xRange / par->N;
    par->dt = par->CFL * par->dx;
    par->Nt = (int)(par->tEND / par->dt);

    par->chanAvg = 0.0;
    par->iterAvg = 0.0;
    par->iterFail = 0;
    par->changeMAX = 0.0;

    if ( (strcmp( par->scheme, "RK" ) == 0) && (par->stage_flag == 0) ) {
        if ( par->M == 2 )
            par->stages = 3;
        else if ( par->M == 3 )
            par->stages = 4;
        else if ( par->M == 4 )
            par->stages = 6;
        else if ( par->M == 5 )
            par->stages = 8;
        else  {
            printf ( "*** ERROR - There is no RK scheme implemented with this order (M=%d)\n", par->M );
            exit ( 1 );
        }
    }

    par->time = 0.0;
    par->count = 0;
    par->N_trouble = 0;
    par->N_neigh = 0;
    par->evalCountF = 0;
    par->evalPrediF = 0;
    par->evalTempF = 0;
    par->evalCountS = 0;
    par->evalPrediS = 0;
    par->evalTempS = 0;
    par->MH2 = par->kappa * par->dx * par->dx;

    par->M_sub = par->M;  // optional

}

// #### #### ####

void Calc_Convergence_Test ( char *file ) {
/*
 * run convergence test from separate program in ../conv/
 * use same parameter file
 */

    char command[255];
    sprintf ( command, "../conv/run_TEST %s", file );
    printf ( "%s\n", command );

    if ( !system ( command ) )
        printf ( "# SUCCESS # - Convergence Test Done.\n" );

}
