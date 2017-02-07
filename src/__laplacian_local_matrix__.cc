#include <octave/oct.h>

DEFUN_DLD (__laplacian_local_matrix__, args, ,
           "Computing the local matrix for the Laplace operator")
{
  octave_value retval;
  
  double hx = args(0).double_value ();
  double hy = args(1).double_value ();
  double D_loc = args(2).double_value ();
  
  Matrix A_loc (4, 4, 0.0);
  
  double diag_coeff = (hx * hx + hy * hy) / (4 * hx * hy);
  
  A_loc.xelem(0, 0) = diag_coeff;
  A_loc.xelem(1, 1) = diag_coeff;
  A_loc.xelem(2, 2) = diag_coeff;
  A_loc.xelem(3, 3) = diag_coeff;
  
  A_loc.xelem(0, 1) = -hy / (2 * hx);
  A_loc.xelem(0, 3) = -hx / (2 * hy);
  A_loc.xelem(1, 2) = -hx / (2 * hy);
  A_loc.xelem(2, 3) = -hy / (2 * hx);
  
  A_loc += A_loc.transpose();
  A_loc *= D_loc;
  
  retval = A_loc;
  return (retval);
};
