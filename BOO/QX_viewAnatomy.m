% This function aims to view the anatomies specified in structures.json
% files.

addpath(genpath('BOO_QL'));
addpath(genpath('CERR2016'));
addpath(genpath('CERRaddins'));
addpath(genpath('utilities'));
addpath(genpath('beamlogs'));

globalFolder = '/data/qifan/dataset_qlyu/UCLAPatients';
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 8;

for i = 1:numPatients
    patientName = ['patient', num2str(i)];
    patExpFolder = fullfile(expFolder, patientName);
    CERRdataFile = fullfile(patExpFolder, [patientName, '.mat']);
    CERR('CERRSLICEVIEWER')
    sliceCallBack_QL('OPENNEWPLANC', CERRdataFile);
    printf(patientName);
end