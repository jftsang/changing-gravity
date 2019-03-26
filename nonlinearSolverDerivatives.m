%% Work out shear rate and inertial number
uzg = zeros(size(ug));
ig = zeros(size(ug));
%{
for tind = 1:nt
  uzg(1, tind) = ( -(3/2)*ug(1, tind)  + 2*ug(2, tind) - (1/2)*ug(3, tind))/dz;
  for j = 2:nz-1
    uzg(j, tind) = ( ug(j+1, tind) - ug(j-1, tind) )/(2*dz);
  end
  uzg(nz, tind) = ( (3/2)*ug(nz, tind) - 2*ug(nz-1, tind) + (1/2)*ug(nz-2, tind))/dz;
  
  ig(:, tind) = uzg(:, tind) ...
                    ./ ( g(ts(tind)) * cos(theta(ts(tind))) * (1-zs') ).^(1/2);
end
%}

for tind = 1:nt
    uzg(1, tind) = (ug(2, tind) - ug(1, tind) ) / dzs(1);
    for j = 2:nz-1
        uzg(j, tind) = ( ug(j+1, tind) - ug(j-1, tind) ) / (dzs(j-1) + dzs(j));
    end
    uzg(nz, tind) = (ug(nz, tind) - ug(nz-1, tind) ) / dzs(nz-1);

    ig(:, tind) = uzg(:, tind) ...
                      ./ ( g(ts(tind)) * cos(theta(ts(tind))) * (1-zs') ).^(1/2);
end
clear tind j;
