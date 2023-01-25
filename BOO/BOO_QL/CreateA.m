function [A,Weights] = CreateA(M, StructureInfo, varargin)
tic
if(~isempty(varargin))
    modFactor = varargin{1};
    CalcA = 1;
    if(length(varargin)==2)
        CalcA = varargin{2};
    end
else
    modFactor = 1;
    CalcA = 1;
end

PTVind = 1;
BODYind = 2;
PTV = StructureInfo(PTVind).Mask;
PTVdilate = imdilate(PTV,ones(3,3,3));

% BODY = StructureInfo(BODYind).Mask;
% M(BODY==0,:) = 0;

StructureName = {StructureInfo.Name};
IsPTV = contains(StructureName,'PTV','IgnoreCase',true);
IsOAR = ~IsPTV;

numVOI = length(StructureInfo);
A_mats = cell(numVOI,1);
maxDose_vecs = cell(numVOI,1);
maxWeights_vecs = cell(numVOI,1);
minDoseTarget_vecs = cell(numVOI,1);
minDoseTargetWeights_vecs = cell(numVOI,1);
OARWeights_vecs = cell(numVOI,1);
numVoxs = zeros(numVOI,1);

IsVOI = true(1,numVOI);
for idx = 1:numVOI
    IsVOI(idx) = IsZeroOrNaN([StructureInfo(idx).maxWeights...
        StructureInfo(idx).minDoseTargetWeights StructureInfo(idx).OARWeights]);
end

PTVsind = find(IsPTV&IsVOI);
OARsind = find(IsOAR&IsVOI);
numPTV = length(PTVsind);
PTV0mask = false(size(PTV));
PTV0mask_dilate = imdilate(PTV0mask,ones(6,6,6));

for idx = 1:numPTV
    StructureInfo(PTVsind(idx)).Mask(PTV0mask_dilate ==1) = 0;
    if idx<numPTV
        PTV0mask = (PTV0mask | StructureInfo(PTVsind(idx)).Mask);
        PTV0mask_dilate = imdilate(PTV0mask,ones(3,3,3));
    end
end

for idx = 1:numVOI
    if(isempty(find(idx==find(IsVOI), 1)))
        continue
    end
    
    if(find(idx==OARsind))
        StructureInfo(idx).Mask(PTVdilate==1) = 0;
    end
    fnd = find(StructureInfo(idx).Mask);
    [I, J, K] = ind2sub(size(PTV),fnd);
    fnd2 = (mod(I,modFactor) == 0 & mod(J,modFactor) == 0 & mod(K,modFactor) == 0);
    fnd = fnd(fnd2);
    
    if(CalcA==1)
        A_mats{idx} = M(fnd,:);
    end
    
    numVox = length(fnd);
    numVoxs(idx) = numVox;
    maxDose_vecs{idx} = StructureInfo(idx).maxDose + zeros(numVox,1); % maximum dose vector
    maxWeights_vecs{idx} = StructureInfo(idx).maxWeights + zeros(numVox,1);
    minDoseTarget_vecs{idx} = StructureInfo(idx).minDoseTarget + zeros(numVox,1);
    minDoseTargetWeights_vecs{idx} = StructureInfo(idx).minDoseTargetWeights + zeros(numVox,1);
    OARWeights_vecs{idx} = StructureInfo(idx).OARWeights + zeros(numVox,1);

end

Weights.maxDose = cat(1,maxDose_vecs{[PTVsind,OARsind]});
Weights.maxWeightsLong = cat(1,maxWeights_vecs{[PTVsind,OARsind]});
Weights.minDoseTarget = cat(1,minDoseTarget_vecs{PTVsind});
Weights.minDoseTargetWeights = cat(1,minDoseTargetWeights_vecs{PTVsind});
Weights.OARWeightsLong = cat(1,OARWeights_vecs{OARsind});


A_ptv = cat(1,A_mats{PTVsind});
A_noPtv = cat(1,A_mats{OARsind});
A = [A_ptv;A_noPtv];

toc

end

function flag = IsZeroOrNaN(x)
flag = ~all(isnan(x)|x==0);
end





