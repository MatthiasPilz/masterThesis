#ifndef ADVECTION_H
#define ADVECTION_H

// prototype declarations

void Advection_Var_Num              ( PARA *par );
void Advection_Init_Data            ( FV_DATA *V, PARA *par );
void Advection_Output_Averages      ( double **V, PARA *par, char *name );
void Advection_Free                 ( FV_DATA *V, PARA *par );
void Advection_Flux_Coeff_i         ( ADER_DATA *A, PARA *par, int i );
void Advection_Flux_Coeff_all       ( ADER_DATA *A, PARA *par );
double Advection_Calc_Speed         ( double **V, PARA *par, int i );
void Advection_Convergence_Output   ( FV_DATA *V, PARA *par, char *name );
void Advection_Phy_Admiss           ( DG_DATA *D, PARA *par );
void Advection_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );

#endif // ADVECTION_H
