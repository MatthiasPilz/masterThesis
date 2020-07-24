// prototype declarations for output subdir

void Output_Parameter               ( PARA *par );
void Setup_Output_Data_Structure    ( PARA *par );
void Output_Polynom_Solution        ( gsl_vector ***W, double *LAM, PARA *par, char *name );
void Output_GSL_Matrix              ( gsl_matrix *A, int sizeR, int sizeC );
void Output_Evolution_Summary       ( PARA *par );
void Output_Current_Status          ( PARA *par );
void Output_Current_Status_Subcell  ( PARA *par );
void Output_Final_Status            ( PARA *par, PARA *par_AW );
void Output_Troubled_Cells          ( int **beta, PARA *par, char *name );
void Output_Troubled_Cells_Conv     ( int **beta, PARA *par, char *name );
void Output_Double_Matrix           ( double **A, int sizeR, int sizeC );
void Output_Troubled_Cells_QS       ( int **beta, PARA *par, char *name );
void Output_Polynom_Solution_on_Lagrange_Points ( gsl_vector ***W, double *LAM, PARA *par, char *name );
