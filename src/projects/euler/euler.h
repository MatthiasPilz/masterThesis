// euler project
void Euler_Var_Num          ( PARA *par );
void Euler_Init_Data		( FV_DATA *V, PARA *par );
void Euler_Output_Averages  ( double **V, PARA *par, char *name );
void Euler_Free             ( FV_DATA *V, PARA *par );
void Euler_Flux_Coeff_i     ( ADER_DATA *A, PARA *par, int i );
void Euler_Flux_Coeff_all   ( ADER_DATA *A, PARA *par );
double Euler_Calc_Speed_av  ( double **V, PARA *par, int i );
void Euler_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );
void Euler_Convergence_Output ( FV_DATA *V, PARA *par, char *name );

//void Wave_Convergence_Output ( FV_DATA *V, PARA *par, char *name );
