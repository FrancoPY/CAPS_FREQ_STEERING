function [ba,scale,phiS,phiR] = calcDeplBa(f1,z,sam,ref,model)

    % Convert the prepared absolute depletion ratio to B/A
    % Supply alpha1/alpha2 in Np/m, c0 in m/s, reference ba, and depth in m
    z = z(:);
    [phiS,normS] = calcDeplPhi(sam.alpha1,sam.alpha2 ./ sam.alpha1,z,model);
    [phiR,normR] = calcDeplPhi(ref.alpha1,ref.alpha2 ./ ref.alpha1,z,model);

    % Cancel the common squared depth before evaluating the surface limit
    scale = sqrt(normR ./ normS) .* (sam.c0 ./ ref.c0).^3 .* (1 + ref.ba / 2);
    ba = 2 * (f1 .* scale - 1);

end

