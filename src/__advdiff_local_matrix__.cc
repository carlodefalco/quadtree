// Copyright (C) 2004-2016  Carlo de Falco
//
// This program is Free Software
//     secs3d - A 3-D Drift--Diffusion Semiconductor Device Simulator
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; If not, see <http://www.gnu.org/licenses/>.
//
//  author: Carlo de Falco

#include <octave/oct.h>

void
bimu_bernoulli (double x, double &bp, double &bn)
{
  const double xlim = 1.0e-2;
  double ax  = fabs (x);

  bp  = 0.0;
  bn  = 0.0;

  //  X=0
  if (x == 0.0)
    {
      bp = 1.0;
      bn = 1.0;
      return;
    }

  // ASYMPTOTICS
  if (ax > 80.0)
    {
      if (x > 0.0)
        {
          bp = 0.0;
          bn = x;
        }
      else
        {
          bp = -x;
          bn = 0.0;
        }
      return;
    }

  // INTERMEDIATE VALUES
  if (ax <= 80 &&  ax > xlim)
    {
      bp = x / (exp (x) - 1.0);
      bn = x + bp;
      return;
    }

  // SMALL VALUES
  if (ax <= xlim &&  ax != 0.0)
    {
      double jj = 1.0;
      double fp = 1.0;
      double fn = 1.0;
      double df = 1.0;
      double segno = 1.0;
      while (fabs (df) > 1.0e-16)
        {
          jj += 1.0;
          segno = -segno;
          df = df * x / jj;
          fp = fp + df;
          fn = fn + segno * df;
        }
      bp = 1 / fp;
      bn = 1 / fn;
      return;
    }

};

DEFUN_DLD (__advdiff_local_matrix__, args, ,
           "Computing the local matrix for the advection-diffusion operator")
{

  octave_value retval;
  
  double hx = args(0).double_value ();
  double hy = args(1).double_value ();
  double alpha_loc = args(2).double_value ();
  Array<double> psi_loc = args(3).array_value ();
  
  double psi12 = psi_loc(1) - psi_loc(0);
  double psi23 = psi_loc(2) - psi_loc(1);
  double psi34 = psi_loc(3) - psi_loc(2);
  double psi41 = psi_loc(0) - psi_loc(3);

  double bp12, bm12, bp23, bm23, bp34, bm34, bp41, bm41;
  bimu_bernoulli (psi12, bp12, bm12);
  bimu_bernoulli (psi23, bp23, bm23);
  bimu_bernoulli (psi34, bp34, bm34);
  bimu_bernoulli (psi41, bp41, bm41);
  
  bp12 = alpha_loc * bp12 * hy / (2 * hx);
  bm12 = alpha_loc * bm12 * hy / (2 * hx);
  bp23 = alpha_loc * bp23 * hx / (2 * hy);
  bm23 = alpha_loc * bm23 * hx / (2 * hy);
  bp34 = alpha_loc * bp34 * hy / (2 * hx);
  bm34 = alpha_loc * bm34 * hy / (2 * hx);
  bp41 = alpha_loc * bp41 * hx / (2 * hy);
  bm41 = alpha_loc * bm41 * hx / (2 * hy);

  Matrix A_loc (4, 4, 0.0);
  
  A_loc.xelem (0, 0) = bm12 + bp41;
  A_loc.xelem (0, 1) = -bp12;
  A_loc.xelem (0, 3) = -bm41;
  
  A_loc.xelem (1, 0) = -bm12;
  A_loc.xelem (1, 1) = bp12 + bm23;
  A_loc.xelem (1, 2) = -bp23;
  
  A_loc.xelem (2, 1) = -bm23;
  A_loc.xelem (2, 2) = bp23 + bm34;
  A_loc.xelem (2, 3) = -bp34;
  
  A_loc.xelem (3, 0) = -bp41;
  A_loc.xelem (3, 2) = -bm34;
  A_loc.xelem (3, 3) = bp34 + bm41;

  retval = A_loc;
  return (retval);
};
