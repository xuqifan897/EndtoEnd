% this code loads and revises the parameters of optimization
% firstly, set the beamweights to 150, and OAR weights to 10, which is
% trail number 1.

if false
globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 8;

%% initial parameter adjustment
trialNO = 2;

% for patientIdx = 1:numPatients
% patientIdx = 5;

list = [3, 5];
for iii = 1:length(list)
    patientIdx = list(iii);
    patientName = ['patient', num2str(patientIdx)];
    patExpFolder = fullfile(expFolder, patientName);
    optFolder = fullfile(patExpFolder, 'optimize');
    paramFile = fullfile(optFolder, 'params0.mat');
    StructInfoFile = fullfile(optFolder, 'StructureInfo0.mat');
    
    load(paramFile);
    load(StructInfoFile);

    params.beamWeight = 1200;
    params.showTrigger = 500;

    nStructures = length(StructureInfo);
    for j = 1:nStructures
        name = StructureInfo(j).Name;
        if strcmpi(name, 'PTV') || strcmpi(name, 'BODY')
            continue
        end
        StructureInfo(j).OARWeights = 10;
    end

    paramsOutFile = fullfile(optFolder, ['params', num2str(trialNO), '.mat']);
    StructureInfoOutFile = fullfile(optFolder, ['StructureInfo', num2str(trialNO), '.mat']);
    save(paramsOutFile, 'params');
    save(StructureInfoOutFile, 'StructureInfo');
end

%%  here is to change the parameters for the optimization. 
% Previously, we compromied the PTV too much.
globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 8;
for i = 1:numPatients
    patientName = ['patient', num2str(i)];
    expPatFolder = fullfile(expFolder, patientName);
    optFolder = fullfile(expPatFolder, 'optimize');

    % find trailNO
    pattern = fullfile(optFolder, 'StructureInfo*.mat');
    files = dir(pattern);
    trailNOs = [];
    for j = 1:length(files)
        digitFlag = isstrprop(files(j).name, 'digit');
        trailNO = str2num(files(j).name(digitFlag));
        trailNOs(end+1) = trailNO;
    end
    maxTrailNO = max(trailNOs);
    newTrailNO = maxTrailNO + 1;

    oldParamsFile = fullfile(optFolder, ['params', num2str(maxTrailNO), '.mat']);
    oldStructureFile = fullfile(optFolder, ['StructureInfo', num2str(maxTrailNO), '.mat']);
    load(oldParamsFile, 'params')
    load(oldStructureFile, 'StructureInfo')

    % PTV weight remain unchanged, change the OAR weight to 2
    % parameters unchanged.
    for j = 1:length(StructureInfo)
        if ~ (strcmpi(StructureInfo(j).Name, 'PTV') || strcmpi(StructureInfo(j).Name, 'BODY'))
            StructureInfo(j).OARWeights = 2;
        end
    end
    newParamsFile = fullfile(optFolder, ['params', num2str(newTrailNO), '.mat']);
    newStructureFile = fullfile(optFolder, ['StructureInfo', num2str(newTrailNO), '.mat']);
    save(newParamsFile, 'params');
    save(newStructureFile, 'StructureInfo');
    patientName
end
end

%% this is to update the parameters to outperform the clinical plans
globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 8;

patientList = [1, 2, 3, 4, 5, 6, 8];
sourceNO = [2, 2, 3, 2, 3, 2, 2, 2];
Names = {{'O_Esgs', 'O_Hrt'}, ...
    {'O_Esgs', 'O_Hrt', 'O_Bwel_Sm', 'O_Bwel_Lg'}, ...
    {'O_Duod', 'O_Esgs', 'O_Hrt'}, ...
    {'O_Hrt', 'O_Stmc'}, ...
    {'O_Cord', 'O_Hrt', 'O_Lung_Tt', 'O_Vssl', 'O_Duod'}, ...
    {'O_Cord', 'O_Duod'}, ...
    {'O_Esgs', 'O_Vessle', 'O_Stmc'}};
Weights = {[5,5], ...
    [5, 5, 5, 5], ...
    [5, 5, 5], ...
    [5, 5], ...
    [5, 5, 5, 3, 5], ...
    [5, 5], ...
    [5, 3, 3]};

% firstly, remove the skin field of patient7
patientName = 'patient7';
optFolder = fullfile(expFolder, patientName, 'optimize');
[paramsFile, StructureFile] = getFiles(optFolder, 2);
load(paramsFile, 'params');
load(StructureFile, 'StructureInfo');
StructureInfo(13).Name
StructureInfo(13) = [];
[paramsOut, StructureOut] = getFiles(optFolder, 3);
save(paramsOut, 'params');
save(StructureOut, 'StructureInfo');

for i = 1:length(patientList)
    idx = patientList(i);
    patientName = ['patient', num2str(idx)];
    names = Names{i};
    weights = Weights{i};
    source = sourceNO(i);
    optFolder = fullfile(expFolder, patientName, 'optimize');
    [ParamsFile, StructureFile] = getFiles(optFolder, source);
    clearvars params StructureInfo
    load(ParamsFile, 'params');
    load(StructureFile, 'StructureInfo');
    StructureInfo = reviseStructures(StructureInfo, names, weights);
    [ParamsOut, StructureOut] = getFiles(optFolder, source+1);
    save(ParamsOut, 'params');
    save(StructureOut, 'StructureInfo');
end

function [paramsFile, StructureFile] = getFiles(optFolder, sourceNO)
    paramsFile = fullfile(optFolder, ['params', num2str(sourceNO), '.mat']);
    StructureFile = fullfile(optFolder, ['StructureInfo', num2str(sourceNO), '.mat']);
end

function StructureInfo = reviseStructures(StructureInfo, names, weights)
    for i = 1:length(StructureInfo)
        for j = 1:length(names)
            if strcmpi(StructureInfo(i).Name, names{j})
                StructureInfo(i).OARWeights = weights(j);
                break
            end
        end
    end
end