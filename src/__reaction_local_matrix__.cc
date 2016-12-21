#include <octave/oct.h>

DEFUN_DLD (__reaction_local_matrix__, args, ,
           "Computing the local matrix for the advection-diffusion operator")
{
  octave_value retval;
  
  double hx = args(0).double_value ();
  double hy = args(1).double_value ();
  double delta_loc = args(2).double_value ();
  
  Matrix A_loc (4, 4, 0.0);
  
  double coeff = delta_loc * hx * hy / 4;
  
  A_loc.xelem (0, 0) = coeff;
  A_loc.xelem (1, 1) = coeff;
  A_loc.xelem (2, 2) = coeff;
  A_loc.xelem (3, 3) = coeff;
  
  retval = A_loc;
  return (retval);
};
