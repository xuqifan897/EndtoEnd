function [xFull,costDFs,costs] = polish_BOO_IMRT_cpu(Beam,A,D,Weights,params)
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

pruneTrigger = 40;
% pruneTrigger = 10000; disp('WARNING!\nWARNING!\nWARNING! PRUNING IS OFF')
t = params.stepSize;
maxIter = params.maxIter;
showTrigger = params.showTrigger;
gamma = params.gamma; % Huber parameter
eta = params.eta;
reductionFactor = .5; % This is a line search parameter.
beamSizes = params.beamSizes;

maxDose = Weights.maxDose;
maxWeightsLong = Weights.maxWeightsLong;
minDoseTarget = Weights.minDoseTarget;
minDoseTargetWeights = Weights.minDoseTargetWeights;
OARWeightsLong = Weights.OARWeightsLong;
numVoxPTV = length(minDoseTarget);

A_ptv = A(1:numVoxPTV,:);
A_noPtv = A(numVoxPTV+1:end,:);

cumsumbeamSizes = cumsum([0;beamSizes]);

beamletselect = zeros(cumsumbeamSizes(end),1);
for i = 1:length(Beam)
    idx = Beam(i);
    beamletselect(cumsumbeamSizes(idx)+1:cumsumbeamSizes(idx+1)) = 1;
end
beamletList = find(beamletselect==1);
A = A(:,beamletList);
D = D(:,beamletList);
ATrans = A';
DTrans = D';
A_ptv = A_ptv(:,beamletList);
A_noPtv = A_noPtv(:,beamletList);

numBeamlets = size(A,2); %numBeamlets is never reduced

    function [grad,cost,costDF] = eval_grad(x)
        
        cost = 0;
        Ax = ATrans'*x;
        
        prox = min(Ax(1:numVoxPTV) - minDoseTarget,0);
        grad = A_ptv'*(minDoseTargetWeights.*(prox)); % mu is still a scalar
        myVal = .5*sum(minDoseTargetWeights.*((prox).^2));
        cost = cost + myVal;
        
        prox = max(Ax - maxDose,0);
        grad = grad + A'*(maxWeightsLong.*(prox));
        myVal = .5*sum(maxWeightsLong.*((prox).^2));
        cost = cost + myVal;
        
        term = Ax(numVoxPTV+1:end);
        grad = grad + A_noPtv'*(OARWeightsLong.*term);
        cost = cost + .5*sum(OARWeightsLong.*(term.^2));
        costDF = cost;
        term = DTrans'*x;
        prox = prox1Norm(term,gamma);
        grad = grad + eta*(D'*(term - prox)/gamma);
        myVal = sum(abs(prox)) + (.5/gamma)*sum((prox - term).^2);
        cost = cost + eta*myVal;
        
    end

    function [cost,costDF] = eval_g(x) % This function computes gx without extra expense of computing grad_gx.
        
        cost = 0;
        Ax = ATrans'*x;
                
        prox = min(Ax(1:numVoxPTV) - minDoseTarget,0);
        myVal = .5*sum(minDoseTargetWeights.*((prox).^2));
        cost = cost + myVal;
        
        prox = max(Ax - maxDose,0);
        myVal = .5*sum(maxWeightsLong.*((prox).^2));
        cost = cost + myVal;
        
        term = Ax(numVoxPTV+1:end);
        cost = cost + .5*sum(OARWeightsLong.*(term.^2));
        costDF = cost;
        term = DTrans'*x;
        prox = prox1Norm(term,gamma);
        myVal = sum(abs(prox)) + (.5/gamma)*sum((prox - term).^2);
        cost = cost + eta*myVal;
        
    end

xkm1 = rand(numBeamlets,1);
vkm1 = xkm1;

disp(['all zero cost is: ',num2str(eval_g(0*xkm1))])

costs = []; costDFs = [];

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
    costs = [costs,cost];
    costDFs = [costDFs,costDF];
        
    if mod(k,showTrigger) == 0
        figure(200)
        semilogy(costs);
    end
    
end

xFull = zeros(cumsumbeamSizes(end),1);
xFull(beamletList) = xkm1;

end
