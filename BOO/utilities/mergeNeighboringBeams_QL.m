function [beamsToKeep,keepIdxs] = mergeNeighboringBeams_QL(ang, activeBeams, numBeamsWeWant)
% This code uses k-means clustering to find 'numBeamsWeWant' beams to use during treatment.

if nargin < 2
    numBeamsWeWant = 20;
end

rng(1)
ang = ang(activeBeams,:);
azimuth = ang(:,1);
elevation = 90-ang(:,2);
r = 1;
[x,y,z] = sph2cart(azimuth,elevation,r);
points = [x,y,z]; 
[clusterIdxs,C] = kmeans(points,numBeamsWeWant,'Replicates',100);
% visualizeBeamClusters(activeBeams,clusterIdxs)

beamsToKeep = [];
keepIdxs = [];
for clusterIdx = 1:numBeamsWeWant
    
    mu = C(clusterIdx,:);
    fnd = find(clusterIdxs == clusterIdx);
    clusterAng = points(fnd,:);
    mnDiff = Inf;
    mnIndx = -1;
    for indx = 1:length(fnd)
        diff = norm(clusterAng(indx,:) - mu);
        if(diff < mnDiff)
            mnDiff = diff;
            mnIndx = indx;
        end
    end    
    keepIdxs = [keepIdxs,fnd(mnIndx)];
    beamsToKeep = [beamsToKeep,activeBeams(fnd(mnIndx))];
                    
end


end