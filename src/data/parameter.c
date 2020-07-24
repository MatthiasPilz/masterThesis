#include "ader.h"
#include "data.h"


//### ### ### ### ### ### ###

void Read_Parameter ( char *name, PARA *par ) {
/*
 * name - filename of parameter file
 * par  - parameter
 *
 * read parameter file and store parameter in par
 * set dependent parameters
 */

    ASSERT ( name );
    sprintf ( par->file, "%s", name );

	char input[256];
	char var[64];
	char wert[64];
	char filename[255];
	
	sprintf ( filename, "../par/%s", name );

	FILE *fp = fopen ( filename, "r" );
    ASSERT ( fp );

// set parameters from file
	while ( fgets ( input, 256, fp ) ) {
		sscanf ( input, "%s\t=\t%s", var, wert );
		if ( strcmp ( var, "M" ) == 0 ) {
			par->M = atoi ( wert );
        }
        else if ( strcmp ( var, "M_sub" ) == 0 ) {
            par->M_sub = atoi ( wert );
        }
		else if ( strcmp ( var, "project" ) == 0 ) {
			sprintf ( par->project, "%s", wert );
		}
        else if ( strcmp ( var, "scheme" ) == 0 ) {
            sprintf ( par->scheme, "%s", wert );
        }
        else if ( strcmp ( var, "stages" ) == 0 ) {
            par->stages = atoi ( wert );
        }
		else if ( strcmp ( var, "N" ) == 0 ) {
            par->N = atoi ( wert );
		}
		else if ( strcmp ( var, "xRange" ) == 0 ) {
			par->xRange = atof ( wert );
		}
        else if ( strcmp ( var, "hack_bound_R" ) == 0 ) {
            par->hack_bound_R = atof ( wert );
        }
		else if ( strcmp ( var, "xMin" ) == 0 ) {
			par->xMin = atof ( wert );
        }
        else if ( strcmp ( var, "CFL_sub" ) == 0 ) {
            par->CFL_sub = atof ( wert );
        }
        else if ( strcmp ( var, "boundary" ) == 0 ) {
            sprintf ( par->boundary, "%s", wert );
            if( (strcmp ( par->boundary, "symmetric" ) == 0) || (strcmp( par->boundary, "hack" ) == 0) ) {
                par->N = par->N / 2 + par->N % 2;
                if(par->xMin + par->xRange < 0) { printf("Grid entirely on the negative x axis is not compatible with symmetry mode!\n"); exit(0); }
                else { par->xRange += par->xMin; par->xMin = 0; }
            }

        }
		else if ( strcmp ( var, "tEND" ) == 0 ) {
			par->tEND = atof ( wert );
		}
        else if ( strcmp ( var, "iterEXTRA" ) == 0 ) {
            par->iterEXTRA = atoi ( wert );
        }
		else if ( strcmp ( var, "eps" ) == 0 ) {
			par->eps = atof ( wert );
		}
		else if ( strcmp ( var, "output" ) == 0 ) {
			par->output = atoi ( wert );
		}
        else if ( strcmp ( var, "Npoly" ) == 0 ) {
            par->Npoly = atoi ( wert );
        }
        else if ( strcmp ( var, "kappa" ) == 0 ) {
            par->kappa = atof ( wert );
        }
        else if ( strcmp ( var, "WENO_R" ) == 0 ) {
            par->WENO_R = atoi ( wert );
        }
        else if ( strcmp ( var, "WENO_eps" ) == 0 ) {
            par->WENO_eps = atof ( wert );
        }
        else if ( strcmp ( var, "WENO_lam" ) == 0 ) {
            par->WENO_lam = atof ( wert );
        }
        else if ( strcmp ( var, "shape" ) == 0 ) {
            sprintf ( par->shape, "%s", wert );
        }
        else if ( strcmp ( var, "dimension" ) == 0 ) {
            par->dimension = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_N_start" ) == 0 ) {
            par->conv_N_start = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_N_end" ) == 0 ) {
            par->conv_N_end = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_M_start" ) == 0 ) {
            par->conv_M_start = atoi ( wert );
        }
        else if ( strcmp ( var, "conv_M_end" ) == 0 ) {
            par->conv_M_end = atoi ( wert );
        }
        else if ( strcmp ( var, "change_M" ) == 0 ) {
            par->change_M = atoi ( wert );
        }
        else if ( strcmp ( var, "change_N" ) == 0 ) {
            par->change_N = atoi ( wert );
        }
        else if ( strcmp ( var, "e_gamma" ) == 0 ) {
            par->e_gamma = atof ( wert );
        }
		else {
			printf ( "*** Can't read parameter %s\n", var );
		}
	}
    fclose ( fp );

// checking for meaningful parameters
    ASSERT ( par->CFL_sub > 0.0 );
    ASSERT ( par->tEND > 0.0 );
    ASSERT ( par->xRange > 0.0 );
    ASSERT ( par->M > 1 && par->M < 7 );
    ASSERT ( par->N > 2*par->M+1 );
    ASSERT ( par->WENO_R > 0 && par->WENO_R < 20 );
    ASSERT ( par->WENO_eps > 0 && par->WENO_eps < 1 );
    ASSERT ( par->WENO_lam > 0 );

/// SET DEPENDENT PARAMETERS
    // grid
	par->size = pow ( par->M+1, 2 );
    par->Mp1 = par->M+1;
    if ( par->M % 2 ) par->n_stencil = 4;
    else par->n_stencil = 3;

    if ( par->M == 1 )
        par->n_stencil = 2;

    par->NSC = 2*par->M + 1;
    par->CFL = par->CFL_sub / ( 2.0*par->M + 1 );
	par->xMax = par->xMin + par->xRange;
	par->dx = (double)par->xRange / par->N;

    // time
    par->dt = par->CFL * par->dx;
	par->Nt = (int) ( par->tEND / par->dt );
	par->time = 0.0;
    par->count = 0;
	
    // iteration
    par->iterMAX = par->Mp1 + par->iterEXTRA;
    par->chanAvg = 0.0;
    par->changeMAX = 0.0;
	par->iterAvg = 0.0;
	par->iterFail = 0;
    par->iterFailCount = 0;

    if ( par->stages == 0 )
        par->stage_flag = 0;
    else
        par->stage_flag = 1;

    // Runge Kutta
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

    /// distinguish between 1D and 3D TOV evolution
    if ( strcmp ( par->shape, "TOV" ) == 0 )
        sprintf ( par->shape, "%s_%dD", par->shape, par->dimension );

    // troubled cells
    par->admiss_num = 0;
    par->admiss_phy = 0;
    par->MH2 = par->kappa * par->dx * par->dx;
    par->N_trouble = 0.0;

    // function evaluation statistic
    par->evalCountF = 0;
    par->evalPrediF = 0;
    par->evalTempF = 0;
    par->evalCountS = 0;
    par->evalPrediS = 0;
    par->evalTempS = 0;
	
    par->var = 0;       // will be set later!

    par->M_sub = par->M;

/// TODO -- not necessary... better: always evolve whole domain and extract special data only in convergence test
//    if ( strcmp( par->project, "grhdc" ) == 0 ) {
//        if ( par->conv_grhdc_R1 <= par->xMax )
//            par->conv_grhdc_N1 = (int)( par->conv_grhdc_R1 / par->dx );
//        else {
//            printf ( "*** ERROR - Radius for grhdc convergence test\n" );
//            exit ( 1 );
//        }
//    }

	Output_Parameter ( par );
}


void Set_ADER_WENO_Parameter ( PARA *D, PARA *A ) {
/*
 * D - parameter on grid
 * A - parameter on subgrid
 *
 * set parameter for subcells according to known grid parameters
 */

    A->M_sub = D->M; // save but don't use!

    sprintf ( A->file, "%s", D->file );
    sprintf ( A->boundary, "%s", D->boundary );
    sprintf ( A->scheme, "%s", D->scheme );
    sprintf ( A->project, "%s", D->project );
    sprintf ( A->shape, "%s", D->shape );

    // grid
    A->M = D->M_sub;
    A->Mp1 = A->M + 1;
    A->size = pow ( A->M+1, 2 );
    if ( A->M % 2 ) A->n_stencil = 4;
    else A->n_stencil = 3;

    if ( A->M == 1 )
        A->n_stencil = 2;

    A->N = D->N * (2*D->M + 1);  // 2M+1 resolution
    A->xRange = D->xRange;
    A->xMin = D->xMin;
    A->xMax = D->xMax;
    A->dx = (double)A->xRange / A->N;
    A->NSC = 2*D->M + 1;
    A->CFL = D->CFL_sub;       // to be able to do one timestep on the subcells
    A->CFL_sub = A->CFL;

    // time
    A->Nt = D->Nt;
    A->tEND = D->tEND;      // !?
    A->dt = A->CFL * A->dx;
    A->time = 0.0;
    A->count = 0;

    // iteration
    A->iterMAX = D->iterMAX;
    A->eps = D->eps;
    A->chanAvg = 0.0;
    A->changeMAX = 0.0;
    A->iterAvg = 0.0;
    A->iterFail = 0;
    A->iterFailCount = 0;

    if ( strcmp( D->scheme, "RK" ) == 0 )
        A->stages = D->stages;

    // troubled cells
    A->admiss_num = 0;
    A->admiss_phy = 0;
    A->N_trouble = 0;
    A->N_neigh = 0;
    A->kappa = D->kappa;
    A->MH2 = A->kappa * A->dx * A->dx;

    A->WENO_eps = D->WENO_eps;
    A->WENO_lam = D->WENO_lam;
    A->WENO_R = D->WENO_R;

    // output
    A->output = 99999;
    A->Npoly  = 0;

    // function evaluation statistic
    A->evalCountF = 0;
    A->evalPrediF = 0;
    A->evalTempF = 0;
    A->evalCountS = 0;
    A->evalPrediS = 0;
    A->evalTempS = 0;

    // convergence
    A->conv_M_end = D->conv_M_end;
    A->conv_M_start = D->conv_M_start;
    A->conv_N_end = D->conv_N_end* (2*D->M + 1);
    A->conv_N_start = D->conv_N_start* (2*D->M + 1);
    A->change_M = D->change_M;
    A->change_N = D->change_N;

    A->e_gamma = D->e_gamma;

    A->var = D->var;
    A->source_flag = D->source_flag;
    A->speed = D->speed;

    A->sym_flags = D->sym_flags;
}

// #### #### ####

void Update_Time_Parameter ( PARA *par, PARA *par_AW ) {

    par->time = par->dt * par->count;
    par_AW->time = par->time;
    par_AW->count = par->count;

}







