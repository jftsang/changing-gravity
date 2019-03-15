%% Bagnoldian solution
%%   function [bagSol, Ithts] 
%%     = nonlinearSolverBagnoldians(mu, g, theta, ts, zs)
function [bagSol, Ithts] = nonlinearSolverBagnoldians(mu, g, theta, ts, zs)
  if (nargin ~= 5)
      print_usage();
  end
  thetas = arrayfun(theta, ts);
  thetas = min(thetas):(0.05*pi/180):max(thetas);

  %% Need to invert mu(I)
  for k = 1:length(thetas)
      if (k == 1)
          Iths(k) = fzero(@(I) mu(I) - tan(thetas(k)), 0);
      else
          Iths(k) = fzero(@(I) mu(I) - tan(thetas(k)), Iths(k-1));
      end
  end
  
  nz = length(zs);
  nt = length(ts);
  [tg, zg] = meshgrid(ts, zs);
  if (length(thetas) == 1)
      Ithts = repmat(Iths, 1, nt);
  else
      Ithts = interp1(thetas, Iths, arrayfun(theta, ts));
  end
  
  uBag = zeros(nz, nt);
  for tind = 1:nt
      uBag(:, tind) = 2/3 * Ithts(tind) * (g(ts(tind)) * cos(theta(ts(tind))))^(1/2) ...
      * ( 1 - (1-zs).^(3/2) );
  end
  
  bagSol = struct('ts', ts, 'zs', zs, 'tg', tg, 'zg', zg, 'ug', uBag);
end
