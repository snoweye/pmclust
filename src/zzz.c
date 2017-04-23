#include <R.h>
#include <R_ext/Rdynload.h>

#include "zzz.h"

static const R_CallMethodDef callMethods[] = {
	{"W_plus_y_k", (DL_FUNC) &W_plus_y_k, 5},
	{"W_plus_y", (DL_FUNC) &W_plus_y, 4},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_pmclust(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_pmclust(). */
