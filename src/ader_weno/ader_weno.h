// prototype declaration for ader_weno subdir

// ### ### ### ader_weno.c
void    Set_Initial_Iteration_Condition_i   ( ADER_DATA *A, gsl_vector ***W, PARA *par, int i );
double  Calc_Change_k                       ( gsl_vector ***Vold, gsl_vector ***Vnew, int i, int k, PARA *par );
void    Final_Step                          ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro );

// ### ### ### finite_volume.c
void Update_Finite_Volume_Scheme    ( FV_DATA *V, PARA *par, int **beta );
void Calculate_Flux_Averages        ( FV_DATA *V, ADER_DATA *A, PARA *par, int **beta, PROJ *pro );
void Calculate_Source_Averages      ( FV_DATA *V, ADER_DATA *A, PARA *par, int **beta, PROJ *pro );

// ### ### ### weno.c
void    Solve_Algebraic_System  ( WENO_DATA *W, PARA *par, int **beta );
void    Calculate_WENO_Coeff    ( WENO_DATA *W, PARA *par, int **beta );
