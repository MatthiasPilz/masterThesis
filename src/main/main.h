// prototypes for main subdir

// ### ### ###
PROJ    Set_Function_Pointer                ( PARA *par );
void    Set_Scheme_Function_Pointer         ( PARA *par );
void    Set_RK_Dense_Output_Function        ( PARA *par );
void    Delete_Data                         ( PARA *par );
void    Setup_Convergence_Data_Structure    ( PARA *par );
void    Delete_old_Iteration_File           ( PARA *par );
void    Copy_Parameter_Data                 ( char *name, char *path );
void    Save_Iteration_Statistic            ( PARA *par, PARA *subcell );

// ### ### ### ader_dg.c
void ADER_DG_Step          ( DG_DATA *D, ADER_DATA *A, PARA *par, PROJ *pro );
void Accept_Unlimited_DG   ( DG_DATA *DG, FV_DATA *V_AW, WENO_DATA *W_AW, ADER_DATA *A_AW, PARA *par, PARA *par_AW, PROJ *pro );

// ### ### ### iteration.c
void Iteration ( ADER_DATA *A, gsl_vector ***W, int **beta, PARA *par, PROJ *pro );

// ### ### ### WRAPPER
int SUBCELL_WENO_TIME_STEP              ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro );
int ALLOC_INIT_ALL_DATA                 ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, DG_DATA *DG, PARA *par, PROJ *pro );
int ALLOC_INIT_ADER_WENO                ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, PARA *par, PROJ *pro );
int SUBCELL_WENO_TIME_STEP_woStatistics ( FV_DATA *V, WENO_DATA *W, ADER_DATA *A, int **beta, PARA *par, PROJ *pro );

// ### ### ### RK
void set_initial_q_RK_i   ( gsl_vector ****q, gsl_vector ***W, PARA *par, int i );
void calc_RK_k            ( gsl_vector ****q, gsl_vector ****k, gsl_matrix *D, double *LAM, gsl_vector **F, gsl_vector **S, PARA *par, int i, int t, PROJ *pro );
void RK_Predictor         ( ADER_DATA *A, gsl_vector ***W, int **beta, PARA *par, PROJ *pro );

