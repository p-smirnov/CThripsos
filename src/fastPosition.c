#include <R.h>
#include <Rinternals.h>
//Copyright https://github.com/joshuaulrich, under GPLv3
SEXP firstCross(SEXP x, SEXP th, SEXP rel, SEXP start, SEXP direction)
{
  int i, int_rel, int_start;
  double *real_x=NULL, real_th;

  if(ncols(x) > 1)
    error("only univariate data allowed");

  /* this currently only works for real x and th arguments
   * support for other types may be added later */
  real_th = asReal(th);
  int_rel = asInteger(rel);
  int_start = asInteger(start)-1;

  /*
   * This is probably written this way to help the compiler in optimizing these internal loops. Therefore, we duplicate for teh reverse option.
   */

  switch(int_rel) {
  case 1:  /* >  */
  real_x = REAL(x);
    for(i=int_start; i<nrows(x); i++)
      if(real_x[i] >  real_th)
        return(ScalarInteger(i+1));
      break;
  case 2:  /* <  */
  real_x = REAL(x);
    for(i=int_start; i<nrows(x); i++)
      if(real_x[i] <  real_th)
        return(ScalarInteger(i+1));
      break;
  case 3:  /* == */
  real_x = REAL(x);
    for(i=int_start; i<nrows(x); i++)
      if(real_x[i] == real_th)
        return(ScalarInteger(i+1));
      break;
  case 4:  /* >= */
  real_x = REAL(x);
    for(i=int_start; i<nrows(x); i++)
      if(real_x[i] >= real_th)
        return(ScalarInteger(i+1));
      break;
  case 5:  /* <= */
  real_x = REAL(x);
    for(i=int_start; i<nrows(x); i++)
      if(real_x[i] <= real_th)
        return(ScalarInteger(i+1));
      break;
  default:
    error("unsupported relationship operator");
  }
  /* return number of observations if relationship is never TRUE */
  return(ScalarInteger(nrows(x)));
}
