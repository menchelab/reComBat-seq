// Adapted from edgeR 4.0.1 (init.cpp)
// Modified by Zhasmina Stoyanova, 2026.
// Changes: removed unused symbol registrations; 
// Original authors: edgeR team (see edgeR package).
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "utils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(compute_apl, 6),

	CALLDEF(fit_levenberg, 11),
	CALLDEF(get_levenberg_start, 6),

	{NULL, NULL, 0}
};

R_CMethodDef all_c_entries[] = {
    {NULL, NULL, 0}
  };

void attribute_visible R_init_reComBatseq(DllInfo *dll) {
	R_registerRoutines(dll, all_c_entries, all_call_entries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

}
