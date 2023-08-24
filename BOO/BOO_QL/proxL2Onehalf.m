function prox = proxL2Onehalf(g0,tau)
% Applies to vector input g0
% prox-operator of f(x) = ||x||_2^q where q = 1/2. Evaluated at vector g0.
% See the paper "Low rank priors for color image segmentation" by Cremers et al
% Section 3.2

    nrmg0 = norm(g0);
    alpha = tau*(nrmg0^(-1.5)); % Possibly Inf
    if alpha > 2*sqrt(6)/9
        tHat = 0;
    else
        root3 = sqrt(3);
        sHat = (2 / root3) * sin((acos(.75 * root3 * alpha) + pi / 2) / 3);
        tHat = sHat^2;      
    end
    
    prox = tHat*g0;                

end