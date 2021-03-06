%% Bagnoldian solution
%%   function [bagSol, Ithts] 
%%     = nonlinearSolverBagnoldians(mu, g, theta, ts, zs)
function [bagSol, Ithts] = nonlinearSolverBagnoldians(mu, g, theta, ts, zs)
  if (nargin ~= 5)
      print_usage();
  end
  thetas = arrayfun(theta, ts);
  thetas = min(thetas):(0.05*pi/180):max(thetas);

  %% Check if the value of tan(theta) falls within the range of mu. 
  %% (We assume that mu is increasing.)
  Imin = 0; muMin = mu(Imin); % TODO try to avoid I=0!
  Imax = 1000; muMax = mu(Imax);

  %% Need to invert mu(I)
  for k = 1:length(thetas)
      % If out of range, then use +-inf.
      if (tan(thetas(k)) > muMax)
          Iths(k) = +inf;
      elseif (tan(thetas(k)) < muMin)
          Iths(k) = -inf;
      else 
          % Otherwise, try solving using fzero, using previous results to help
          % if possible.
          if (k == 1)
              Iths(k) = fzero(@(I) mu(I) - tan(thetas(k)), 0);
          else
              if ( isinf(Iths(k-1)) )
                  Iths(k) = fzero(@(I) mu(I) - tan(thetas(k)), 0);
              else
                  Iths(k) = fzero(@(I) mu(I) - tan(thetas(k)), Iths(k-1));
              end
          end
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

  dz = zs(2) - zs(1);
  qs = integrateg(uBag, dz);
  chis = integrateg(uBag.^2, dz) ./ qs.^2;

  bagSol = struct('ts', ts, 'dt', ts(2)-ts(1), 'tmax', ts(end), ...
            'zs', zs, 'dz', dz, ...
            'gs', arrayfun(g, ts), 'thetas', arrayfun(theta, ts), ...
            'Ithts', Ithts, ...
            'tg', tg, 'zg', zg, 'ug', uBag, ...
            'qs', qs, 'chis', chis);
end
