## Copyright (C) 2017 Carlo de Falco

clear all
close all

KT  = 1e-7;
Lx  = 0;
Rx  = 45;
mu  = 1;
nu  = 2;
g   = 30;
PM  = 30;


G0  = @(p) (200/pi) * atan (4*(PM - p));
G   = @(p) 1 * ifelse (G0(p) >= 0, G0(p), 0*p);

pressure = @(m, n) ((g + 1)/g) .* (n + m) .^ g;

N     = 225;
T     = 1;
x     = linspace (Lx, Rx, N+1) .';

msh = msh2m_quadtree (x, x);
msh = bim2c_quadtree_mesh_properties (msh, [], []);
Nnodes = columns(msh.p);
Nelements = columns(msh.t);
x = msh.p(1, :).';
y = msh.p(2, :).';

m = .1 * exp (- .5 * (x.^2 + y.^2)) + 1e-9;
n  = .8 * exp (-5e-7 * (x.^2 + y.^2)) + 1e-9;

dnodes = []; %msh2m_nodes_on_sides(msh, 1:4);
inodes = setdiff (1:Nnodes, dnodes);

zel = zeros (Nelements, 1);
zno = zeros (Nnodes, 1);
oel = ones (Nelements, 1);
ono = ones (Nnodes, 1);

mass   = bim2a_quadtree_reaction (msh, oel, ono);

hm = @(u) 4 ./  sum (1 ./ u(msh.t (1:4, :)), 1);
diffm  = @(m, n) hm (mu * g * ((g + 1)/g) .* m .* (n + m) .^ (g-1));
diffn  = @(m, n) hm (nu * g * ((g + 1)/g) .* n .* (n + m) .^ (g-1));
stiffm = @(m, n) bim2a_quadtree_laplacian (msh, diffm (m,n)); 
stiffn = @(m, n) bim2a_quadtree_laplacian (msh, diffn (m,n)); 

function [res, jac] = fun (u, mold, nold, inodes, dnodes,
                           stiffm, stiffn, mass, pressure, G, dt)

  mest = mold;
  nest = nold;
  mest(inodes) = u(1:numel(inodes));
  nest(inodes) = u(numel(inodes)+1:end);

  Sm  = stiffm (mest, nest);
  Sn  = stiffn (mest, nest);
  GG  = mass .* sparse (diag (G (pressure (mest, nest))));
    
  Amm = (mass + dt * (Sm - GG));
  Amn = dt * Sm;  
  bm  = mass * mold;
    
  Ann = (mass + dt * Sn);
  Anm = dt * Sn;
  bn  = mass * nold;

  A = [Amm(inodes, inodes), Amn(inodes, inodes);
       Anm(inodes, inodes), Ann(inodes, inodes)];

  b = [(bm(inodes) - Amm(inodes,dnodes) * mold(dnodes, :) -
        Amn(inodes,dnodes) * nold(dnodes, :));
       (bn(inodes) - Anm(inodes,dnodes) * mold(dnodes, :) -
        Ann(inodes,dnodes) * nold(dnodes, :))];

  res = A*u-b;
  jac = A;
  
endfunction

nt    = 100;
tsave = linspace (0, T, nt);
ii    = 2;
t     = 0;
tstore= t;
dt    = 1e-5;

mvold = mold = m;
nvold = nold = n;

for its = 2 : nt
  
  while (t < tsave(its))

    if (t + dt > tsave(its))
      dt = tsave(its) - t;      
    endif
    
    printf ("t = %g, dt = %g\n", t, dt)

    mvold = mold;
    nvold = nold;
    
    mold = m;
    nold = n;

    uguess = [mold;nold];
    if (ii > 3)
      mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      if (! all ([mguess;nguess] >= 0))
        printf ("negative guess\n")
      endif

      uguess = max (0, [mguess;nguess]);
      
    endif
    
    while (true)

      uold = [mold;nold]; 
      u = fsolve (@(u) ...
                   fun (u, mold, nold, inodes, dnodes,
                        stiffm, stiffn, mass, pressure, G, dt),
                  uguess, optimset ("Jacobian", "on"));

      all_nzero = all (u >= 0);
      incr = norm ((uguess - u)./(u + 1), inf);
      small_err = incr < 1e-4;

      if (all_nzero && small_err)
        break
      else
        dt *= .5;
        printf (["reducing time step : t = %g, ", ...
                 "dt = %g, all_nzero = %d, small_err = %d\n"],
                t, dt, all_nzero, small_err)
        uguess = [mold;nold];      
        if (ii > 3)
          mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          if (! all ([mguess;nguess] >= 0))
            printf ("negative guess\n")
          endif
          uguess = max (0, [mguess;nguess]);
          
        endif

      endif
    
    endwhile

  
    t += dt;
    tstore(ii) = t;
    if (incr == 0)
      dt *= 2;
    else
      dtfact = min (sqrt(.38) * sqrt (1e-4 / incr), 2);
      dt = min (dtfact * dt, .1);
    endif
  
    m = mold;
    m(inodes) = u(1:numel(inodes));

    n = nold;
    n(inodes) = u(numel(inodes)+1:end);

    assert (all (n >= 0))
    assert (all (m >= 0))
  
    ## figure (1)
    ## plot (x, m(:, [ii-1:ii]), x, n(:, [ii-1:ii]))
    ## title (sprintf ("t = %g", t))
    ## figure (2)
    ## plot (x, pressure (m(:, ii), n(:, ii)))

    ## figure (3)
    ## plot ((x(2:end) + x(1:end-1))/2, zz (x, m(:, ii)),
    ##       (x(2:end) + x(1:end-1))/2, zz (x, n(:, ii)))
    
    ## drawnow
    ii += 1;

  endwhile

  filename = sprintf ("test_tumor_%4.4d", its);
  delete ([filename ".vtu"])
  fpl_vtk_write_field_quadmesh (filename, msh, {m, "m"; n, "n"}, {}, 1);
  printf ("saving step : t = %g its = %d\n", t, its)
  
endfor



