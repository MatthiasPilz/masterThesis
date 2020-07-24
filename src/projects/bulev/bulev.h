//### ### wave project:
void BuLev_Var_Num           ( PARA *par );
void BuLev_Init_Data         ( FV_DATA *V, PARA *par );
void BuLev_Free              ( FV_DATA *V, PARA *par );
void BuLev_Flux_Coeff_i      ( ADER_DATA *A, PARA *par, int i );
void BuLev_Flux_Coeff_all    ( ADER_DATA *A, PARA *par );
double BuLev_Calc_Speed      ( double **V, PARA *par, int i );

void BuLev_Convergence_Output ( FV_DATA *V, PARA *par, char *name );
