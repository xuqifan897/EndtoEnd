function [prox, nrmnew] = proxL2Onehalf_gpu(g0,tau)
% Applies to matrix input g0; evaluate proximal operator for each column of
% matrix g0
% prox-operator of f(x) = ||x||_2^q where q = 1/2. Evaluated at each column of
% matrix g0
% See the paper "Low rank priors for color image segmentation" by Cremers et al
% Section 3.2
global root3 root6

nrm2 = sum(g0.^2,1);
nrm234 = nrm2(:).^(-3/4);
tnrm234 = tau.*nrm234;
alpha = gather(tnrm234); % Possibly Inf
sHat = (2 / root3) * sin((acos(.75 * root3 * alpha) + pi / 2) / 3);
tHat = sHat.^2;
tHat(alpha > 2*root6/9) = 0;
tHat = gpuArray(tHat);
prox = (tHat'.*g0);

nrm2newbuff = sum(prox.^2,1);
nrmnew = sqrt(nrm2newbuff);
nrmnew = gather(nrmnew);

end