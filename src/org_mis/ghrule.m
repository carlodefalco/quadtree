## Gauss-Hermite quadrature rule.
function [z, p] = ghrule (n)
  if (n <= 1)
    p = sqrt (pi);
    z = 0;
    return
  end

  jac = zeros (n);
  v = sqrt ((1:n) / 2);
  jac = jac + diag (v(1:n-1), 1) + diag (v(1:n-1), -1);
  [p, z] = eig (jac);
  norm2 = sqrt (diag (p' * p));              # Weight normalization.
  p = (sqrt (pi) * p(1,:)' .^ 2) ./ norm2;   # sqrt(pi) = beta0;
  z = diag (z);
endfunction
