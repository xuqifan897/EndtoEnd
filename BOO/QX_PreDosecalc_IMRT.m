% This file runs on shenggpu2. Dose calculation preprocessing
addpath(genpath('BOO_QL'));
addpath(genpath('CERR2016'));
addpath(genpath('CERRaddins'));
addpath(genpath('utilities'));
addpath(genpath('beamlogs'));

globalFolder = '/data/qifan/dataset_qlyu/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
beamlogfile = 'spine_T1-6_beamlog.mat';
PrescriptionDose = 20;
numPatients = 8;

if ~ isfolder(expFolder)
    mkdir(expFolder)
end

for i = 2:numPatients
    patientName = ['patient', num2str(i)];
    patDataFolder = fullfile(dataFolder, patientName);
    patExpFolder = fullfile(expFolder, patientName);
    if ~ isfolder(patExpFolder)
        mkdir(patExpFolder);
    end
    
    %% load to CERR
    DicomPath = fullfile(patDataFolder, 'CTAlignResize');
    OutputFileName = fullfile(patExpFolder, [patientName, '.mat']);
    CERRImportDicom_QL(DicomPath, OutputFileName);
    CERR('CERRSLICEVIEWER')
    sliceCallBack_QL('OPENNEWPLANC', OutputFileName);
    
    %% generate beam_list file
    load('4pi_angles.mat');
    load(beamlogfile);
    
    % for each patient, we choose half of the beams by random
    select = probabilitySampling(size(theta_VarianIEC, 1), 0.5);
    
    MLCangle = 0; % MLC angle is set to zero;
    Gantry = theta_VarianIEC(select,1);
    Couch = theta_VarianIEC(select,2);
    MLCangles = MLCangle*ones(length(Gantry),1);
    Angles = [Gantry Couch MLCangles];
    
    beamlistFile = fullfile(patExpFolder, 'beamlist.txt');
    fileID = fopen(beamlistFile, 'w');
    for ii = 1:size(Angles, 1)
        fprintf(fileID,'%6.4f %6.4f %6.4f \n',Angles(ii,:));
    end
    fclose(fileID);
    
    %% generate structures.json file
    baseFileNames = dir([DicomPath, '/*.dcm']);
    count2 = 0;
    for ii = 1:length(baseFileNames)
        FileName = baseFileNames(ii).name;
        fullFileName = fullfile(DicomPath, FileName);
        info = dicominfo(fullFileName);
        if(strcmp(info.Modality, 'RTSTRUCT'))
            RTstructurefiles = FileName;
            if(count2 > 1)
                error('Multiple RT structure files!');
            end
            count2 = count2 + 1;
        end
    end
    RTstructInfo = dicominfo(fullfile(DicomPath, RTstructurefiles));
    allstructs = fieldnames(RTstructInfo.StructureSetROISequence);
    for ii = 1:length(allstructs)
        structures{ii} = RTstructInfo.StructureSetROISequence.(...
            allstructs{ii}).ROIName;
    end
    jsonFileName = fullfile(patExpFolder, 'structures.json');
    SaveStructureFileOnly(structures, jsonFileName, PrescriptionDose);
end

function select = probabilitySampling(size, probability)
    select = logical(zeros(size, 1));
    numSelected = floor(size * probability);
    select(1:numSelected) = 1;
    select = select(randperm(size));
end