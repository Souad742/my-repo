/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab,incx,incy;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double alpha,beta;
  double *RHS, *EX_SOL, *X;
  double *AB;
  double temp, relres;

  NRHS=1;
  nbpoints=8;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  //lab=kv+kl+ku+1;
  lab=kv+kl+ku;
  incx=1;
  incy=2;
  alpha=1.0;
  beta=0.0;

  AB = (double *) malloc(sizeof(double)*lab*la);

  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 1; 


  if (row == 0){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    
      //info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
      


  } 
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

  // info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
   //cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, EX_SOL, incx, beta, RHS, incy);
  
  }    

    const double b[3]={3.0,1.0,1.0};
 const double x[3]={3.0,3.0,2.0};
 double x_s[3]={0.0,0.0,0.0};
 const double Bb[3][3] = {{ 2 , -1 , 0}, { -1 , 2 , -1 }, { 0, -1 , 2}};
  
  
  int i;
  int j;

   for(i=0;i<3;i++){
   for(j=0;j<3;j++)
     { x_s[i]=x_s[i]+(alpha*Bb[i][j]*x[j])+(beta*x_s[i]);}
     printf("%f\n",x_s[i]);  
       
    }
      cblas_dgbmv(CblasColMajor, CblasNoTrans, 3, lab, kl, ku, alpha, Bb, 4, x, incx, beta, x_s, incy);
        for(i=0;i<3;i++){
       printf("%f\n",x_s[i]); } 





  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");
  //write_vec(RHS, &la, "RHS.dat");
 
  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nThe relative residual error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
}
