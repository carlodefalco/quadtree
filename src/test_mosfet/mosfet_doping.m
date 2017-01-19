function [Nd, Na] = mosfet_doping (msh)

  N_plus  = 1e25;
  N_minus = 1e24;
  P_plus  = 1e25;
  P_minus = 1e24;

  [x, y] = deal (msh.p(1, :).', msh.p(2, :).');

  xl = min (x);
  xr = max (x);
  L  = xr - xl;
  
  yb = min (y);
  yt = max (y);
  H  = yt - yb;
  
  Na  = P_minus * ones (size (x));
  Na += P_plus  * gaussian (y, -H/2, H/10) .* ...
        (((x>L/3) & (x<2*L/3)) + ...
         gaussian (x, L/3, L/10) .* (x <= L/3) + ...
         gaussian (x, 2*L/3, L/10) .* (x >= 2*L/3));

  Nd  = zeros (size (x));
  Nd += N_plus  * gaussian (y, 0, H/4) .* ...
        ((x<L/4) + ...
         (x>3*L/4) + ...
         gaussian (x, L/4, L/10) .* (x >= L/4) + ...
         gaussian (x, 3*L/4, L/10) .* (x <= 3*L/4));
  
endfunction

function y = gaussian (x, c, s)
  y = exp (- (x-c) .^2 ./ (2 * s^2));
endfunction
