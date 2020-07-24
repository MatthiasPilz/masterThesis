//### ### wave project:
void Wave_Var_Num           ( PARA *par );
void Wave_Init_Data         ( FV_DATA *V, PARA *par );
void Wave_Output_Averages   ( double **V, PARA *par, char *name );
void Wave_Free              ( FV_DATA *V, PARA *par );
void Wave_Flux_Coeff_i      ( ADER_DATA *A, PARA *par, int i );
void Wave_Flux_Coeff_all    ( ADER_DATA *A, PARA *par );
void Wave_Source_Coeff_i    ( ADER_DATA *A, PARA *par, int i );
void Wave_Source_Coeff_all  ( ADER_DATA *A, PARA *par );
double Wave_Calc_Speed      ( double **V, PARA *par, int i );

void Wave_Convergence_Output ( FV_DATA *V, PARA *par, char *name );

void Wave_S_RK ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM );
void Wave_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );
