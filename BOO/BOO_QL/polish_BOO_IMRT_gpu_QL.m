function [xFull,costDFs,costs] = polish_BOO_IMRT_gpu_QL(Beam,A_cpu,D_cpu,Weights,params)
% This version attempts to increase t less often.

% minimize sum_{b=1}^{numBeams} beamWeights(b) * || x_b ||_2
% .5*mu || (A_0 x - minDose_target)_- ||_2^2
% + sum_{i=0}^numOars .5*alphai || (A_i x - d_vecs{i})_+ ||_2^2
% + sum_{i=1}^numOars .5*betai || A_i x ||^2 + eta || Dx ||_1^gamma
% subject to x >= 0.

% activeBeams contains the NEW indices of the active beams.
% So you could never have 1162 be an active beam, because 1162 is not one
% of the NEW indices.
% (The new indices range from 1 to numValidBeams.  Of the 1162 beams, only
% some are valid in the sense of not causing collisions.

global root3 root6
root3 = sqrt(3);
root6 = sqrt(6);

t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;
gamma = params.gamma; % Huber parameter
eta = params.eta;
reductionFactor = .5; % This is a line search parameter.
beamSizes = params.beamSizes;
numBeams = length(beamSizes);

maxDose_cpu = Weights.maxDose;
maxWeightsLong_cpu = Weights.maxWeightsLong;
minDoseTarget_cpu = Weights.minDoseTarget;
minDoseTargetWeights_cpu = Weights.minDoseTargetWeights;
OARWeightsLong_cpu = Weights.OARWeightsLong;
numVoxPTV = length(minDoseTarget_cpu);

cumsumbeamSizes = cumsum([0;beamSizes]);
beamletselect = zeros(cumsumbeamSizes(end),1);
for i = 1:length(Beam)
    idx = Beam(i);
    beamletselect(cumsumbeamSizes(idx)+1:cumsumbeamSizes(idx+1)) = 1;
end
beamletList = find(beamletselect==1);
A_cpu = A_cpu(:,beamletList);
D_cpu = D_cpu(:,beamletList);

A = gpuArray(A_cpu);
ATrans = A';
D = gpuArray(D_cpu);
DTrans = D';

maxWeightsLong = gpuArray(maxWeightsLong_cpu);
maxDose = gpuArray(maxDose_cpu);
minDoseTarget = gpuArray(minDoseTarget_cpu);
minDoseTargetWeights = gpuArray(minDoseTargetWeights_cpu);
OARWeightsLong = gpuArray(OARWeightsLong_cpu);

numBeamlets = size(A,2); %numBeamlets is never reduced

    function [grad,cost] = eval_grad(x)
        
        Ax = A*x;
        prox1 = min(Ax(1:numVoxPTV) - minDoseTarget,0);
        prox2 = max(Ax - maxDose,0);
        term3 = Ax(numVoxPTV+1:end);
        term4 = D*x;
        prox4 = prox1Norm(term4,gamma);
        
        grad = ATrans*([minDoseTargetWeights.*(prox1);OARWeightsLong.*term3] ...
            + maxWeightsLong.*(prox2)) + eta*(DTrans*(term4 - prox4)/gamma);
        cost = .5*sum(minDoseTargetWeights.*((prox1).^2)) + ...
            .5*sum(maxWeightsLong.*((prox2).^2)) + .5*sum(OARWeightsLong.*(term3.^2))...
            + eta*(sum(abs(prox4)) + (.5/gamma)*sum((prox4 - term4).^2));
        
    end

    function [cost,costDF] = eval_g(x) % This function computes gx without extra expense of computing grad_gx.
        
        Ax = A*x;
        prox1 = min(Ax(1:numVoxPTV) - minDoseTarget,0);
        prox2 = max(Ax - maxDose,0);
        term3 = Ax(numVoxPTV+1:end);
        term4 = D*x;
        prox4 = prox1Norm(term4,gamma);
        
        costDF = .5*sum(minDoseTargetWeights.*((prox1).^2)) +...
            .5*sum(maxWeightsLong.*((prox2).^2)) + .5*sum(OARWeightsLong.*(term3.^2));        
        cost = costDF + eta*(sum(abs(prox4)) + (.5/gamma)*sum((prox4 - term4).^2));
        
    end

xkm1 = gpuArray(rand(numBeamlets,1));
vkm1 = xkm1;
disp(['all zero cost is: ',num2str(eval_g(0*xkm1))])
costs = zeros(maxIter,1);
costDFs = zeros(maxIter,1);

for k = 1:maxIter
    
    if (k <= 50 || mod(k,5) == 0)
        t = t/reductionFactor; % Attempt to increase t.
    end
    accept_t = 0;
    while accept_t == 0
        if k > 1
            a = tkm1;
            b = t*theta_km1^2;
            c = -t*theta_km1^2;
            
            theta = (-b + sqrt(b^2 - 4*a*c))/(2*a);
            y = (1 - theta)*xkm1 + theta*vkm1;
        else
            theta = 1;
            y = xkm1;
        end
        [gradAty,gy] = eval_grad(y);
        
        in = y - t*gradAty;
        x = max(in,0);
        
        [gx,costDF] = eval_g(x);
        lhs = gx;
        rhs = gy + gradAty'*(x - y) + (.5/t)*sum((x-y).^2);
        if lhs <= rhs
            accept_t = 1;
        else
            t = reductionFactor*t;
        end
        
    end
    
    v = xkm1 + (1/theta)*(x - xkm1);
    
    theta_km1 = theta;
    tkm1 = t;
    xkm1 = x;
    vkm1 = v;
    
    cost = gx;    
    costs(k) = gather(cost);
    costDFs(k) = gather(costDF);
        
    if mod(k,showTrigger) == 0
        figure(200)
        semilogy(costs);
    end
    
end

xFull = zeros(cumsumbeamSizes(end),1,'gpuArray');
xFull(beamletList) = xkm1;

end
