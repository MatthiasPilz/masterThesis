//### ### grhdc project:
void grhdc_Var_Num                     ( PARA *par );
void grhdc_Init_Data                   ( FV_DATA *V, PARA *par );
void grhdc_Output_Averages             ( double **V, PARA *par, char *name );
void grhdc_Output_Averages_single_time ( double **V, PARA *par, char *name );
void grhdc_Free                        ( FV_DATA *V, PARA *par );
void grhdc_Flux_Coeff_i                ( ADER_DATA *A, PARA *par, int i );
void grhdc_Flux_Coeff_all              ( ADER_DATA *A, PARA *par );
void grhdc_Source_Coeff_i_1D           ( ADER_DATA *A, PARA *par, int i );
void grhdc_Source_Coeff_i_3D           ( ADER_DATA *A, PARA *par, int i );
void grhdc_Source_Coeff_all            ( ADER_DATA *A, PARA *par );
double grhdc_Calc_Speed                ( double **V, PARA *par, int i );
void grhdc_Convergence_Output          ( FV_DATA *V, PARA *par, char *name );
void grhdc_F_RK                        ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );
void grhdc_S_RK_1D                     ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM );
void grhdc_S_RK_3D                     ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM );
void grhdc_update_primitives_double_all(double **vec, PARA *par);
void grhdc_update_primitives_all       (gsl_vector ***vec, int length, PARA *par);
void grhdc_update_primitives_i         (gsl_vector ***vec, int length, PARA *par, int i);
void grhdc_update_primitives_i_t       (gsl_vector ****vec, int length, PARA *par, int i, int t);
void grhdc_Phy_Admiss                   ( DG_DATA *D, PARA *par );

//### ### TOV hack for the inside of the star
void Hack_Temp_TOV_Flux                 ( ADER_DATA *A, FV_DATA *V, ADER_DATA *A_AW, FV_DATA *V_AW, PARA *par, int temp_hack_N );
void Hack_Set_DG_Parameter              (PARA *par);
void Hack_Copy_TOV_Flux                 ( ADER_DATA *Afrom, ADER_DATA *Ato, FV_DATA *FVfrom, FV_DATA *FVto, PARA *par );
