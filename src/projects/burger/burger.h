//### ### wave project:
void Burger_Var_Num                     ( PARA *par );
void Burger_Init_Data                   ( FV_DATA *V, PARA *par );
void Burger_Output_Averages             ( double **V, PARA *par, char *name );
void Burger_Output_Averages_single_time ( double **V, PARA *par, char *name );
void Burger_Free                        ( FV_DATA *V, PARA *par );
void Burger_Flux_Coeff_i                ( ADER_DATA *A, PARA *par, int i );
void Burger_Flux_Coeff_all              ( ADER_DATA *A, PARA *par );
void Burger_Source_Coeff_i              ( ADER_DATA *A, PARA *par, int i );
void Burger_Source_Coeff_all            ( ADER_DATA *A, PARA *par );
double Burger_Calc_Speed                ( double **V, PARA *par, int i );
void Burger_Convergence_Output          ( FV_DATA *V, PARA *par, char *name );
void Burger_F_RK                        ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );
void Burger_S_RK                        ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM );
