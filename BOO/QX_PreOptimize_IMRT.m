% this code runs on CPU server, which requires large storage
% firstly, take a try. We previously divided the whole coefficient matrix
% into several partitions. Let's see if it works correctly

addpath(genpath('BOO_QL'));
addpath(genpath('CERR2016'));
addpath(genpath('CERRaddins'));
addpath(genpath('utilities'));
addpath(genpath('beamlogs'));

globalFolder = '/media/raid0/qifan/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 1;
nPartitions = 6;
thresh = 1e-06;

for idx = 1:numPatients
    patientName = ['patient', num2str(idx)];
    patExpFolder = fullfile(expFolder, patientName);
    optFolder = fullfile(patExpFolder, 'optimize');
    if ~ isfolder(optFolder)
        mkdir(optFolder)
    end
    maskfile = fullfile(patExpFolder, 'Dose_Coefficients.mask');
    DoseCoeFolder = fullfile(patExpFolder, 'Dose_Coefficients');
    for j = 1:nPartitions
        parFile = fullfile(DoseCoeFolder, ['partition', num2str(j), '.h5']);
        [M,dose_data,masks]=BuildDoseMatrix(parFile, maskfile, thresh);
        outFile = fullfile(DoseCoeFolder, ['partition', num2str(j), '.mat']);
        save(fullfile(DoseCoeFolder, ['partition', num2str(j), '.mat']),'M','dose_data','masks','-v7.3');
        break
    end
end