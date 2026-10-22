function [phi,normPhi] = calcDeplPhi(a,r,z,model)

    % Supply a in Np/m, r = alpha2/alpha1, and z in m
    % Return phi and phi/(a*z)^2, whose surface limit is 1/8

    % Expand the attenuation and depth over compatible property arrays
    eta = a .* z + zeros(size(r));
    r = r + zeros(size(eta));
    normPhi = zeros(size(eta));
    switch model
        case 'quadratic'
            ratio = -expm1(-2 * eta) ./ eta;
            ratio(eta == 0) = 2;
            normPhi = ratio.^2 / 32;
        case 'general'
            % Evaluate the small-depth series without subtracting similar losses
            small = max(abs(2 * eta),abs(r .* eta)) <= 0.5;
            u = 2 * eta(small);
            v = r(small) .* eta(small);
            poly = ones(size(u));
            vPow = ones(size(u));
            total = 0.5 * ones(size(u));
            fact = 2;
            for k = 1:18
                vPow = vPow .* v;
                poly = u .* poly + vPow;
                fact = fact * (k + 2);
                total = total + (-1)^k * poly / fact;
            end
            normPhi(small) = total / 4;

            % Evaluate the removable attenuation-ratio limit continuously
            near = ~small & abs(r - 2) < 0.01;
            x = eta(near);
            delta = r(near) - 2;
            total = zeros(size(x));
            for k = 0:8
                total = total + (-delta / 2).^k .* gammainc(2 * x,k + 2);
            end
            normPhi(near) = total ./ (16 * x.^2);

            % Evaluate the closed form away from both removable limits
            other = ~small & ~near;
            x = eta(other);
            ratio = r(other);
            phi = (0.5 * -expm1(-2 * x) - -expm1(-ratio .* x) ./ ratio) ./ (4 * (ratio - 2));
            quad = ratio == 4;
            phi(quad) = (-expm1(-2 * x(quad))).^2 / 32;
            normPhi(other) = phi ./ x.^2;
        otherwise
            error('Choose general or quadratic for the depth model.')
    end
    phi = eta.^2 .* normPhi;

end
