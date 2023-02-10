% this code loads and revises the parameters of optimization
% firstly, set the beamweights to 150, and OAR weights to 10, which is
% trail number 1.

globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 8;
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