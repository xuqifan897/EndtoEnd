function params = InitIMRTparams_noInfo(M,dose_data)
%% Beamlet Information
beamletnumTable0 = double([dose_data.beam_metadata.N_beamlets])';
numBeam = length(beamletnumTable0);

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

beamWeightsInit = findBeamWeights(M,beamletnumTable0);
beamWeightsInit = beamWeightsInit/max(beamWeightsInit(:));
beamWeightsInit(beamWeightsInit<0.1) = 0.1;

%% IMRT params
beamWeight = 1; % Adjust this number to control number of active beams
eta = .01; % Smooth parameter
gamma = 1;  % Huber parameter
numBeamsWeWant = 20;
maxIter = 3000;
ChangeWeightsTrigger = 1000;
showTrigger = 10;
tInit = 1e-5;

params = struct('beamWeight',beamWeight,'gamma',gamma,'eta',eta,'numBeamsWeWant',numBeamsWeWant,...
    'stepSize',tInit,'maxIter',maxIter,'showTrigger',showTrigger,...
    'ChangeWeightsTrigger',ChangeWeightsTrigger,'beamWeightsInit',beamWeightsInit,...
    'beamSizes',beamletnumTable0,...
    'BeamletLog0',BeamletLog0);

