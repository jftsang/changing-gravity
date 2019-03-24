function ress = integrateg(vg, dz)
    ress = zeros(1, size(vg, 2));
    for tind = 1:size(vg, 2)
        ress(tind) = integrate(vg(:, tind), dz);
    end
end
