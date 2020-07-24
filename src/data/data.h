// prototype declarations for data subdir

// ### ### ### data.c
void    Read_Gsl_Matrix                 ( gsl_matrix *A, const int sizeR, const int sizeC, char *fullPath );
void    Read_Double_Matrix              ( double **A, const int sizeR, const int sizeC, char *fullPath );
void    Read_Double_Vector              ( double *A, const int size, char *fullPath );
void    Read_GSL_Vector                 ( gsl_vector *A, const int size, char *fullPath );
void    FV_2_WENO_DATA                  ( FV_DATA *V, WENO_DATA *W, PARA *par, int **beta );
void    Set_Beta_OK                     ( int **beta, PARA *par );
void    Set_Beta_Bad                    ( int **beta, PARA *par );
void    Set_Beta_Prev                   ( int **beta, int **beta_prev, PARA *par );
void    Set_Beta_Sub                    ( int **beta, int **beta_sub, PARA *par );
void    Adapt_Troubled_Cells            ( int **beta, PARA *par );
void    Adapt_Troubled_Cells_extended   ( int **beta, PARA *par );
void    Count_Troubled_Cells            ( int *beta, PARA *par );

// ### ### ### parameter.c
void    Read_Parameter              ( char *name, PARA *par );
void    Set_ADER_WENO_Parameter     ( PARA *D, PARA *A );
void    Update_Time_Parameter       ( PARA *par, PARA *par_AW );

// ### ### ### alloc.c
void    Allocate_FV_DATA_general    ( FV_DATA *V, PARA *par );
void    Allocate_WENO_DATA_general  ( WENO_DATA *W, PARA *par );
void    Allocate_ADER_DATA_general  ( ADER_DATA *A, PARA *par );
void    Allocate_DG_DATA_general    ( DG_DATA *DG, PARA *par );

// ### ### ### init.c
void    Init_FV_DATA_general        ( FV_DATA *V, PARA *par );
void    Init_WENO_DATA_general      ( WENO_DATA *W, PARA *par );
void    Init_ADER_DATA_general      ( ADER_DATA *A, PARA *par );
void    Init_DG_DATA_general        ( FV_DATA *V, WENO_DATA *W, DG_DATA *DG, PARA *par, PROJ *pro );

// ### ### ### free.c
void    Free_FV_DATA                ( FV_DATA *V, PARA *par );
void    Free_WENO_DATA              ( WENO_DATA *W, PARA *par );
void    Free_ADER_DATA              ( ADER_DATA *A, PARA *par );
void    Free_DG_DATA                ( DG_DATA *D, PARA *par );

// ### ### ### convergence.c
void    Set_Convergence_Parameter   ( PARA *par );
void    Calc_Convergence_Test       ( char *file );


// ### ### ### timer.c
void    timer_print     ( PARA *par );
int     timer_stop      ( char *name );
int     timer_start     ( char *name );
void    timer_enable    ( int on );
