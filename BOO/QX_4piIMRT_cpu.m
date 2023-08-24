% this code runs optimization
addpath(genpath('BOO_QL'), '-end');
addpath(genpath('CERR2016'), '-end');
addpath(genpath('CERRaddins'), '-end');
addpath(genpath('utilities'), '-end');
addpath(genpath('beamlogs'), '-end');

% users have to specify the parameters in the folloing block
% patientIdx and trailNO can be specified from the command line
if ~ exist('patientIdx', 'var')
    patientIdx = 1;
end
if ~ exist('trailNO', 'var')
    trailNO = 1;
end
globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');

patientName = ['patient', num2str(patientIdx)];
patExpFolder = fullfile(expFolder, patientName);
optFolder = fullfile(patExpFolder, 'optimize');

% patDataFile = fullfile(patExpFolder, [patientName, '.mat']);
% CERR('CERRSLICEVIEWER')
% sliceCallBack_QL('OPENNEWPLANC', patDataFile);

paramsFile = fullfile(optFolder, ['params', num2str(trailNO), '.mat']);
StructureInfoFile = fullfile(optFolder, ['StructureInfo', num2str(trailNO), '.mat']);
matrixFile = fullfile(patExpFolder, ['patient', num2str(patientIdx), '_M.mat']);
matrixOrgFile = fullfile(patExpFolder, ['patient', num2str(patientIdx), '_M_original.mat']);
load(paramsFile, 'params');
load(StructureInfoFile, 'StructureInfo');
load(matrixFile, 'M');
load(matrixOrgFile, 'dose_data', 'masks');

% load the matrix
DS = 1;  % DS = 1 for no downsampling; DS > 1 for downsampling with a factor of DS
[A, Weights] = CreateA(M, StructureInfo, DS);
ATrans = A';

[Dx, Dy] = CreateDxDyFMO(params.BeamletLog0);
D = [Dx; Dy];

% random number
seed = 2;
rng(seed)

%% beam selection
tic
params.showTrigger = 500;
[xFista, costsFista, activeBeams, activeNorms, topN] = BOO_IMRT_L2OneHalf_cpu_QL(A, ATrans, D, Weights, params);
timeBeamSelect = toc;
figure; semilogy(costsFista);
%% debug starts here
BOOresult = struct('patientName',patientName,...
    'params',params,'StructureInfo',StructureInfo,'xFista',xFista,...
    'activeBeams',activeBeams,'activeNorms',activeNorms,...
    'costsFista',costsFista,'timeBeamSelect',timeBeamSelect);

BOOresultFile = fullfile(optFolder, ['BOOresult', num2str(trailNO), '.mat']);
save(BOOresultFile, 'BOOresult');

% show selected beams
finalBeams = activeBeams;
finalBeamsVarianIEC = params.beamVarianIEC(finalBeams, :);
gantryVarianIEC = finalBeamsVarianIEC(:, 1);
couchVarianIEC = finalBeamsVarianIEC(:, 2);

PTV = StructureInfo(1).Mask;
BODY = StructureInfo(2).Mask;
draw_beammask_QL(params.beamfpangles(finalBeams, :), BODY, PTV);

%% Polish step
paramsPolish = params;
paramsPolish.maxIter = 500;
tic
[xPolish, costsDF_polish, costs_polish] = polish_BOO_IMRT_cpu(finalBeams,A,D,Weights,paramsPolish);
timePolish = toc;
figure; semilogy(costsDF_polish);

%% Visualize and save results
dose = M * xPolish;
dose = reshape(dose, size(PTV));
dose(BODY == 0 & PTV == 0) = 0;
polishResult = struct('patientName', patientName, 'dose', dose, 'finalBeams', finalBeams, ...
    'xPolish', xPolish, 'timePolish', timePolish, 'StructureInfo', StructureInfo, ...
    'gantryVarianIEC', gantryVarianIEC, 'couchVarianIEC', couchVarianIEC, ...
    'paramsPolish', paramsPolish);
polishResultFile = fullfile(optFolder, ['polishResult', num2str(trailNO), '.mat']);
save(polishResultFile, 'polishResult');

selected_angles = struct('gantryVarianIEC', gantryVarianIEC,'couchVarianIEC', couchVarianIEC);
T = struct2table(selected_angles);
tableFile = fullfile(optFolder, ['selectedAngles', num2str(trailNO), '.csv']);
writetable(T, tableFile);