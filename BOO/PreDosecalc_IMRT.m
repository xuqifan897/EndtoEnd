clear; clc; close all
patientName = 'spine02';
patFolder = fullfile('D:\datatest\OtherProjects\4piSpine\',patientName);
PresriptionDose = 20;
beamlogfile = 'spine_T1-6_beamlog.mat';

%% Load to CERR
DicomPath = fullfile(patFolder,'dicomdata');
OutputFileName = fullfile(patFolder,[patientName '.mat']);
CERRImportDicom_QL(DicomPath,OutputFileName);
CERR('CERRSLICEVIEWER')
sliceCallBack_QL('OPENNEWPLANC', OutputFileName);

%% generate beam_list file
load('4pi_angles.mat');
load(beamlogfile);

MLCangle = 0; % MLC angle is set to zero;
Gantry = theta_VarianIEC(beamlog_iso==1,1);
Couch = theta_VarianIEC(beamlog_iso==1,2);
MLCangles = MLCangle*ones(length(Gantry),1);
Angles = [Gantry Couch MLCangles];

beamlistFile = fullfile(patFolder,'beamlist.txt');
fileID = fopen(beamlistFile,'w');
for ii = 1:size(Angles,1)
    fprintf(fileID,'%6.4f %6.4f %6.4f \n',Angles(ii,:));
end
fclose(fileID);

%% generate structures.json file
baseFileNames = dir([DicomPath '\*.dcm']);
count2 = 0;
for ii = 1:length(baseFileNames)
    FileName = baseFileNames(ii).name;
    fullFileName = fullfile(DicomPath, FileName);
    info = dicominfo(fullFileName);
    if(strcmp(info.Modality,'RTSTRUCT'))
        RTstructfiles = FileName;
        if(count2>1)
            error('Multiple RT structure files!')
        end
        count2 = count2 + 1;
    end
end
RTstructInfo = dicominfo(fullfile(DicomPath, RTstructfiles));
allstructs = fieldnames(RTstructInfo.StructureSetROISequence);
for ii = 1:length(allstructs)
    structures{ii} = RTstructInfo.StructureSetROISequence.(allstructs{ii}).ROIName;
end

jsonFileName = fullfile(patFolder,'structures.json');
SaveStructureFileOnly(structures, jsonFileName, PresriptionDose);

