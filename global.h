/**********************************************************************
                        Global Variables
**********************************************************************/
char *optstring = "i:o:g:r:s:v:L:a:D:d:I:S:EFCHpmP:";

#if QUICK_FS
double   temps[] = { -1.e-5, -0.075, -0.20, -0.50, -1.00, -2.50, -5, -10 };
#endif

int flag;

global_param_struct global_param;
veg_lib_struct *veg_lib;
option_struct options;
#if LINK_DEBUG
debug_struct debug;
#endif
Error_struct Error;
param_set_struct param_set;
