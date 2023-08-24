function [StructureInfo, params] = InitIMRTparams(M,dose_data,masks,pdose)

StructureNum = length(masks);
for i=1:StructureNum
    StructureName{i} = masks{i}.name;
end

PTVind = 1;
BODYind = 2;
PTV = masks{PTVind}.mask;

masks{PTVind}.name = 'PTV';
masks{BODYind}.name = 'BODY';

[mask_ringstruct, mask_skin] = CreateRingStructureSkin(masks{PTVind}.mask, masks{BODYind}.mask);
RingStructureInd = length(masks) + 1;
masks(RingStructureInd) = cell({struct('name','RingStructure','mask',mask_ringstruct)});
SkinInd = RingStructureInd+1;
masks(SkinInd) = cell({struct('name','Skin','mask',mask_skin)});

IsPTV = contains(StructureName,'PTV','IgnoreCase',true);
if(nnz(IsPTV)>1)
    NotPTV = ~IsPTV;
    IsPTV(PTVind) = 0;
else
    IsPTV(PTVind) = 1;
    NotPTV = ~IsPTV;
end
PTVsind = find(IsPTV);
OARsind = find(NotPTV);

for i = 1:length(OARsind)
    if(OARsind(i)~=BODYind)
        masks{OARsind(i)}.mask = (masks{OARsind(i)}.mask & ~PTV);
    end
end

for i=1:length(masks)
    StructureName{i} = masks{i}.name;
    StructuremaxDose{i} = min(pdose)*0.9;
    StructureminTargetDose{i} = NaN;
    StructureIdealDose{i} = 0;
    StructuremaxWeights{i} = 0;
    StructureWeights{i} = 1;
    StructureminDoseTargetWeights{i} = NaN;
    StructureMask{i} = masks{i}.mask;
    StructureVoxelNum{i} = nnz(masks{i}.mask);
end

StructuremaxDose(PTVsind) = num2cell(pdose);
StructureminTargetDose(PTVsind) = num2cell(pdose);
StructureIdealDose(PTVsind) = num2cell(pdose);
StructuremaxDose(RingStructureInd) = num2cell(min(pdose)*0.9);
StructuremaxWeights(RingStructureInd) = {0};
StructuremaxWeights(PTVsind) = {100};
StructuremaxWeights(BODYind) = {0};
StructuremaxWeights(SkinInd) = {0};
StructureWeights(PTVsind) = {NaN};
StructureWeights(PTVind) = {NaN};
StructureWeights(BODYind) = {0};

StructureminDoseTargetWeights(PTVsind) = {100};

StructureInfo = struct('Name',StructureName,...
    'maxWeights',StructuremaxWeights,'maxDose',StructuremaxDose,...
    'minDoseTargetWeights',StructureminDoseTargetWeights,'minDoseTarget',StructureminTargetDose,...
    'OARWeights',StructureWeights,...
    'IdealDose',StructureIdealDose,'Mask',StructureMask,...
    'VoxelNum',StructureVoxelNum);


%% Beamlet Information
load('4pi_angles','angles','theta_VarianIEC');
beamletnumTable0 = double([dose_data.beam_metadata.N_beamlets])';
numBeam = length(beamletnumTable0);

beamVarianIEC = zeros(numBeam,2);
beamfpangles = zeros(numBeam,2);
for kk = 1:numBeam
    VarianIECgantry = dose_data.beam_metadata(kk).beam_specs.gantry_rot_rad/2/pi*360;
    VarianIECcouch = dose_data.beam_metadata(kk).beam_specs.couch_rot_rad/2/pi*360;  
    beamVarianIEC(kk,:) = [VarianIECgantry,VarianIECcouch];
    [~, minI] = min(sum((beamVarianIEC(kk,:) - theta_VarianIEC).^2,2));
    beamfpangles(kk,:) = angles(minI,:);
end

FmapDim = dose_data.beam_metadata(1).beam_specs.fmap_dims;
[xInd,yInd]=ind2sub(FmapDim,dose_data.column_labels(:,2) + 1);
xlength = max(xInd) - min(xInd) + 1;
ylength = max(yInd) - min(yInd) + 1;

FullFmapDim_new = [xlength,ylength,numBeam];
I_new = xInd - min(xInd) + 1;
J_new = yInd - min(yInd) + 1;
K_new = dose_data.column_labels(:,1) + 1;
Subs=sub2ind(FullFmapDim_new,I_new,J_new,K_new);

BeamletLog0 = zeros(FullFmapDim_new);
BeamletLog0(Subs) = 1;

beamWeightsInit = findBeamWeights(M,beamletnumTable0,PTV);
beamWeightsInit = beamWeightsInit/max(beamWeightsInit(:));
beamWeightsInit(beamWeightsInit<0.1) = 0.1;

%% IMRT params
beamWeight = 50; % Adjust this number to control number of active beams
eta = .1; % Smooth parameter
gamma = 1;  % Huber parameter
numBeamsWeWant = 20;
maxIter = 8000;
ChangeWeightsTrigger = 1000;
showTrigger = 10;
tInit = 1e-5;

params = struct('beamWeight',beamWeight,'gamma',gamma,'eta',eta,'numBeamsWeWant',numBeamsWeWant,...
    'stepSize',tInit,'maxIter',maxIter,'showTrigger',showTrigger,...
    'ChangeWeightsTrigger',ChangeWeightsTrigger,'beamWeightsInit',beamWeightsInit,...
    'beamSizes',beamletnumTable0,...
    'BeamletLog0',BeamletLog0,'beamVarianIEC',beamVarianIEC,'beamfpangles',beamfpangles);


%% Sanity check
cumsumbeamsize = [0; cumsum(params.beamSizes)];
id = find(params.beamVarianIEC(:,2)==0);
nb = 7;
beamsid = randperm(numel(id),nb);
beams = id(beamsid);
for bid = 1:nb
    beam = beams(bid);
    x = ones(size(M,2),1);
    x(cumsumbeamsize(beam)+1:cumsumbeamsize(beam+1)) = 1;
    dose = M*x; dose = reshape(dose, size(PTV));
    PTVdose = dose(PTV==1);
    mask = (dose<median(PTVdose)/3*2)&PTV;
    if(find(mask))
        figure;imshow3D(mask)
        error('FOV might be too small for the target!!!') %% Comment this line if you believe the calculated beamlets are enough to cover the PTV is large enough
    end
end
