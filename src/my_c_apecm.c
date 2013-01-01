#include <R.h>
#include <Rinternals.h>


/*
   This function adds y to W by the i_k-th column and returns W.
   All W, and y should be in double.
*/
SEXP W_plus_y_k(SEXP R_W, SEXP R_y, SEXP R_nrow, SEXP R_ncol, SEXP R_i_k){
	SEXP R_Z;
	double *C_Z, *C_W, *C_y;
	int *C_i_k, *C_nrow, i_k;
	int i;

	C_W = REAL(R_W);
	C_y = REAL(R_y);
	C_i_k = INTEGER(R_i_k);
	C_nrow = INTEGER(R_nrow);

        PROTECT(R_Z = allocVector(REALSXP, *C_nrow));
        C_Z = REAL(R_Z);

	i_k = *C_i_k - 1;
	C_W = C_W + i_k * *C_nrow;
	C_y = C_y + i_k;

	for(i = 0; i < *C_nrow; i++){
		*C_Z = *C_W + *C_y;
		C_Z++;
		C_W++;
	}

	UNPROTECT(1);
	return(R_Z);
} /* End of W_plus_y_k(). */

