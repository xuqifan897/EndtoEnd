function beamWeights = findBeamWeights(M,beamletcount,PTV)

numCandBeams = length(beamletcount);
numBeamlets = size(M,2);
beamWeights = zeros(numCandBeams,1);

idxStop = 0;
for beamIdx = 1:numCandBeams
    
    nBeamlets_thisBeam = beamletcount(beamIdx);
    idxStart  = idxStop + 1;
    idxStop = idxStop + nBeamlets_thisBeam;
    
    x = zeros(numBeamlets,1);
    x(idxStart:idxStop) = 1;
    
    tic
%     Mx = M*x;
    Mx = M(:,idxStart:idxStop)*ones(nBeamlets_thisBeam,1);
    timeMx = toc
    
    beamWeights(beamIdx) = sqrt(mean(Mx(PTV==1))/sqrt(nBeamlets_thisBeam));
   
%     doseAll = reshape(Mx,size(PTV));
% 
%     currentTime = fix(clock);
%     timeLabel = [num2str(currentTime(4)),'_',num2str(currentTime(5)),'_',num2str(currentTime(6))];
%     title(num2str(currentTime))
%     addDoseToGui_dvo(doseAll,['fista_',timeLabel])
    
end







end