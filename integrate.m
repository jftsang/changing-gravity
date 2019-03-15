function res = integrate(vs, dz)
    res = ( 0.5 * vs(1) + sum(vs(2:end-1)) + 0.5 * vs(end) ) * dz;
end
