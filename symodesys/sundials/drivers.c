#include <stdlib.h> // malloc, free
#include "drivers.h"
#include "symodesys_util.h"

#include <stdio.h> // printf

enum {
  STATUS_FOUT = 1000,
  STATUS_Y,
  STATUS_ABSTOL_,
  STATUS_CVODE_MEM,
  STATUS_DKY_OUT,
};

/* 
TODO: 
  make use of: 
      CVodeGetNumRhsEvals
      CVodeGetNumJacEvals <-- cvdls
      
*/


/* Scalar reltol and vector abstol */
int integrate_fixed_step(
   double t, double t1, double * y0, int n_steps,
   double h_init, double h_max, double abstol,
   double reltol, void * params, DIM_T dim,
   int nderiv, double * tout, double * Yout,
   int step_type_idx, int mode, int mu, int ml)
{
  /* 
     Adapted from cvRoberts_dns.c
     step_type_idx:
        1 => CV_ADAMS
	2 => CV_BDF
     mode:
        0 => DENSE_MODE
	1 => BANDED_MODE
  */

  double ti;
  double dt = (t1-t)/(double)(n_steps-1);
  int status = 0;
  N_Vector y = NULL;
  N_Vector abstol_ = NULL;
  void *cvode_mem = NULL;
  N_Vector * dky_out = NULL;


  y = N_VMake_Serial(dim, y0);
  if (y == NULL){
    status = STATUS_Y;
    goto exit_y;
  }


  abstol_ = N_VNew_Serial(dim);
  for (int i=0; i<dim; ++i) NV_Ith_S(abstol_, i) = abstol; // TODO: let user input vector
  if (abstol_ == NULL){
    status = STATUS_ABSTOL_;
    goto exit_abstol_;
  }


  if (step_type_idx == 1) step_type_idx = CV_ADAMS;
  if (step_type_idx == 2) step_type_idx = CV_BDF;  
  // For now we skip CV_FUNCTIONAL only use CV_NEWTON
  cvode_mem = CVodeCreate(step_type_idx, CV_NEWTON); 
  if (cvode_mem == NULL){
    status = STATUS_CVODE_MEM;
    goto exit_cvode_mem;
  }

  // Allocate memory for temporary arrays holding derivatives
  dky_out = malloc(nderiv*sizeof(*dky_out));
  if (dky_out == NULL){
    status = STATUS_DKY_OUT;
    goto exit_dky_out;
  }
  for (int i=0; i<nderiv; ++i){
    dky_out[i] = N_VNew_Serial(dim);
    for (int j=0; j<dim; ++j) NV_Ith_S(dky_out[i], j) = 0.0;
  }
  

  status = CVodeInit(cvode_mem, func, t, y);
  if (status < 0) goto exit_runtime;

  status = CVodeSVtolerances(cvode_mem, reltol, abstol_);
  if (status < 0) goto exit_runtime;

  /* Call CVodeRootInit to specify the root function g with 2 components */
  /* flag = CVodeRootInit(cvode_mem, 2, g); */
  /* if (check_flag(&flag, "CVodeRootInit", 1)) return(1); */

  /* Call CVDense/CVLapackDense to specify the dense linear solver */
  switch(mode){
  case(DENSE_MODE):
    status = OUR_DENSE(cvode_mem, dim);
    if (status < 0) goto exit_runtime;
    /* /\* Set the Jacobian routine to Jac (user-supplied) *\/ */
    status = CVDlsSetDenseJacFn(cvode_mem, dense_jac); 
    break;
  case(BANDED_MODE):
    status = OUR_BAND(cvode_mem, dim, mu, ml);
    if (status < 0) goto exit_runtime;
    status = CVDlsSetBandJacFn(cvode_mem, band_jac); 
    break;
  }
  if (status < 0) goto exit_runtime;

  status = CVodeSetUserData(cvode_mem, params);
  if (status < 0) goto exit_runtime;

  if (h_init > 0.0) CVodeSetInitStep(cvode_mem, h_init);
  if (h_max > 0.0) CVodeSetMaxStep(cvode_mem, h_max);

  /* Run integration */
  for (int i = 0; i < n_steps; ++i){
    printf("i: %i\n", i); fflush(stdout);
    tout[i] = t;
    if (i > 0) 
      for (int k = 0; k < nderiv; ++k)
	CVodeGetDky(cvode_mem, t, k+1, dky_out[k]);
    for (size_t j = 0; j < dim; ++j){
      printf("j: %i\n", j); fflush(stdout);
      Yout[i*dim*(nderiv+1)+j*(nderiv+1)+0] = NV_Ith_S(y, j);
      for (int k=0; k < nderiv; ++k){
	printf("k: %i\n", k); fflush(stdout);
    	Yout[i*dim*(nderiv+1)+j*(nderiv+1)+k+1] = NV_Ith_S(dky_out[k], j);
      }
    }
    /* Macro-step loop */
    ti = t + dt;
    status = CVode(cvode_mem, ti, y, &t, CV_NORMAL); // CV_ONE_STEP: single internal step
    printf("A!\n"); fflush(stdout);

    /* If checking for root it may be printed here */
    if (status != CV_SUCCESS)
      break;
  }

  // Error handling
 exit_runtime:
  CVodeFree(&cvode_mem);  
  for (int i=0; i<nderiv; ++i) N_VDestroy_Serial(dky_out[i]);
 exit_dky_out:
  free(dky_out);
 exit_cvode_mem: 
  N_VDestroy_Serial(abstol_);
 exit_abstol_:
  N_VDestroy_Serial(y);
 exit_y:
  return status;
}

int
integrate_fixed_step_print(double t, double t1, double * y, int n_steps,
			   double h_init, double h_max, double abstol,
			   double reltol, void *params, DIM_T dim, int nderiv,
			   int step_type_idx, int mode, int mu, int ml)
{
  int status;
  double * tout;
  double * Yout;

  tout = malloc(sizeof(double)*n_steps);
  Yout = malloc(sizeof(double)*n_steps*dim*(nderiv+1));
  status = integrate_fixed_step(t, t1, y, n_steps, h_init, h_max, abstol,
				reltol, params, dim, nderiv, tout, Yout, 
				step_type_idx, mode, mu, ml);
  if (status != CV_SUCCESS)
    {
      printf ("Error, return value=%d\n", status);
      goto cleanup;
    }

  for (int i = 0; i < n_steps; ++i)
    {
      print_state(tout[i], dim, nderiv, i, Yout);
    }

 cleanup:
  free(tout);
  free(Yout);
  return status;
}
