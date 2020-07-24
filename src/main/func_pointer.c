#include "ader.h"
#include "main.h"

PROJ Set_Function_Pointer ( PARA *par ) {
	
    PROJ pro;

    pro.f_var_num = NULL;
    pro.f_init_data = NULL;
    pro.f_flux_coeff_i = NULL;
    pro.f_flux_coeff_all = NULL;
    pro.f_source_coeff_i = NULL;
    pro.f_source_coeff_all = NULL;
    pro.f_calc_speed = NULL;
    pro.f_post_evolution = NULL;
    pro.f_convergence_output = NULL;
    pro.f_phy_admiss = NULL;
    pro.f_calc_F_RK = NULL;
    pro.f_calc_S_RK = NULL;
    pro.f_output_averages = NULL;


	if ( strcmp ( par->project, "wave" ) == 0 ) {
        pro.f_var_num               = Wave_Var_Num;
        pro.f_init_data             = Wave_Init_Data;
        pro.f_flux_coeff_i          = Wave_Flux_Coeff_i;
        pro.f_flux_coeff_all        = Wave_Flux_Coeff_all;
        pro.f_source_coeff_i        = Wave_Source_Coeff_i;
        pro.f_source_coeff_all      = Wave_Source_Coeff_all;
        pro.f_calc_speed            = Wave_Calc_Speed;
        pro.f_convergence_output    = Wave_Convergence_Output;
        pro.f_calc_F_RK             = Wave_F_RK;
        pro.f_calc_S_RK             = Wave_S_RK;
        pro.f_output_averages       = Wave_Output_Averages;

        par->source_flag = 1;
        par->speed = sqrt(2.0)/2.0;

    }
    else if ( strcmp ( par->project, "advection" ) == 0 ) {
        pro.f_var_num               = Advection_Var_Num;
        pro.f_init_data             = Advection_Init_Data;
        pro.f_flux_coeff_i          = Advection_Flux_Coeff_i;
        pro.f_flux_coeff_all        = Advection_Flux_Coeff_all;
        pro.f_calc_speed            = Advection_Calc_Speed;
        pro.f_convergence_output    = Advection_Convergence_Output;
        pro.f_phy_admiss            = Advection_Phy_Admiss;
        pro.f_calc_F_RK             = Advection_F_RK;
        pro.f_output_averages       = Advection_Output_Averages;

        par->source_flag = 0;
        par->speed = -sqrt(2.0)/2.0;
    }
    else if ( strcmp ( par->project, "euler" ) == 0 ) {
        pro.f_var_num               = Euler_Var_Num;
        pro.f_init_data             = Euler_Init_Data;
        pro.f_flux_coeff_i          = Euler_Flux_Coeff_i;
        pro.f_flux_coeff_all        = Euler_Flux_Coeff_all;
        pro.f_calc_speed            = Euler_Calc_Speed_av;
        pro.f_convergence_output    = Euler_Convergence_Output;
        pro.f_calc_F_RK             = Euler_F_RK;
        pro.f_output_averages       = Euler_Output_Averages;

        par->source_flag = 0;

    }
    else if ( strcmp ( par->project, "burger" ) == 0 ) {
        pro.f_var_num               = Burger_Var_Num;
        pro.f_init_data             = Burger_Init_Data;
        pro.f_flux_coeff_i          = Burger_Flux_Coeff_i;
        pro.f_flux_coeff_all        = Burger_Flux_Coeff_all;
        pro.f_source_coeff_i        = Burger_Source_Coeff_i;
        pro.f_source_coeff_all      = Burger_Source_Coeff_all;
        pro.f_calc_speed            = Burger_Calc_Speed;
        pro.f_convergence_output    = Burger_Convergence_Output;
        pro.f_calc_F_RK             = Burger_F_RK;
        pro.f_calc_S_RK             = Burger_S_RK;
        pro.f_output_averages       = Burger_Output_Averages;

        par->source_flag = 0;   // optional

    }
    else if ( strcmp ( par->project, "grhdc" ) == 0 ) {
        pro.f_var_num               = grhdc_Var_Num;
        pro.f_init_data             = grhdc_Init_Data;
        pro.f_flux_coeff_i          = grhdc_Flux_Coeff_i;
        pro.f_flux_coeff_all        = grhdc_Flux_Coeff_all;
        pro.f_source_coeff_all      = grhdc_Source_Coeff_all;
        pro.f_calc_speed            = grhdc_Calc_Speed;
        pro.f_convergence_output    = grhdc_Convergence_Output;
        pro.f_phy_admiss            = grhdc_Phy_Admiss;
        pro.f_calc_F_RK             = grhdc_F_RK;
        pro.f_output_averages       = grhdc_Output_Averages;

        if ( strcmp ( par->shape, "TOV_3D" ) == 0 && (strcmp ( par->boundary, "symmetric" ) == 0 || (strcmp ( par->boundary, "hack" ) == 0 )) ) {
            pro.f_calc_S_RK             = grhdc_S_RK_3D;
            pro.f_source_coeff_i        = grhdc_Source_Coeff_i_3D;
            printf ( "*** ATTENTION - set 3D source terms!\n" );
        }
        else {
            pro.f_calc_S_RK             = grhdc_S_RK_1D;
            pro.f_source_coeff_i        = grhdc_Source_Coeff_i_1D;
            printf ( "*** ATTENTION - set 1D source terms!\n" );
        }

        par->source_flag = 1;   // optional

    }
	else {
		printf ( "*** ERROR - function pointer - undefined project!\n*** ABORT!\n" );
		exit ( 1 );
	}
	
	return pro;
	
}

// ### ### ### ### ### ### ###

void Set_Scheme_Function_Pointer ( PARA *par ) {

    if ( strcmp ( par->scheme, "dumbser" ) == 0 )
        f_predictor_step = Iteration;
    else if ( strcmp ( par->scheme, "RK" ) == 0 )
        f_predictor_step = RK_Predictor;
    else {
        printf ( "**** ERROR - undefined predictor scheme!\n*** ABORT!\n" );
        exit ( 1 );
    }
}








