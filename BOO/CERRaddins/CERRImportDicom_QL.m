function CERRImportDicom_QL(DicomPath,bz2FileName)
% Import all dicom files in DicomPath and its subfolders
% save planC in bz2FileName
 
global stateS
if ~isfield(stateS,'initDicomFlag')
    flag = init_ML_DICOM;
    if ~flag
        return;
    end
elseif isfield(stateS,'initDicomFlag') && ~stateS.initDicomFlag
    return;
end
 
tic;
% Read options file
pathStr = getCERRPath;
optName = [pathStr 'CERROptions.m'];
optS = opts4Exe(optName);

hWaitbar = waitbar(0,'Scanning Directory Please wait...');
CERRStatusString('Scanning DICOM directory');
 
dcmdirS = []; 
patientNum = 1;

[filesInCurDir,dirsInCurDir] = rdir(DicomPath);

if isfield(optS,'importDICOMsubDirs') && strcmpi(optS.importDICOMsubDirs,'yes') && ~isempty(dirsInCurDir)    
    
    for i = 1:length(dirsInCurDir)
        %     patient = scandir_mldcm(fullfile(path, dirs(i).name), hWaitbar, i);
        patient = scandir_mldcm(dirsInCurDir(i).fullpath, hWaitbar, i);
        if ~isempty(patient)
            for j = 1:length(patient.PATIENT)
                dcmdirS.(['patient_' num2str(patientNum)]) = patient.PATIENT(j);
                patientNum = patientNum + 1;
            end
        end
    end
    
else
    
    filesV = dir(DicomPath);
    disp(DicomPath);
    dirs = filesV([filesV.isdir]);
    dirs(2) = [];
    dirs(1).name = '';
    
    for i = 1:length(dirs)
        patient = scandir_mldcm(fullfile(DicomPath, dirs(i).name), hWaitbar, i);
        if ~isempty(patient)
            for j = 1:length(patient.PATIENT)
                dcmdirS.(['patient_' num2str(patientNum)]) = patient.PATIENT(j);
                patientNum = patientNum + 1;
            end
        end
    end
    
end


if isempty(dcmdirS)
    close(hWaitbar);
    msgbox('There is no dicom data!','Application Info','warn');
    return;
end
 
close(hWaitbar);
 
selected = 'all';
patNameC = fieldnames(dcmdirS);
if isempty(selected)
    return
elseif strcmpi(selected,'all')
    combinedDcmdirS = struct('STUDY',dcmdirS.(patNameC{1}).STUDY,'info',dcmdirS.(patNameC{1}).info);
    count = 0;
    for studyCount = 1:length(combinedDcmdirS.STUDY)
        for seriesCount = 1:length(combinedDcmdirS.STUDY(studyCount).SERIES)
            count = count + 1;
            newCombinedDcmdirS.STUDY.SERIES(count) = combinedDcmdirS.STUDY(studyCount).SERIES(seriesCount);
        end
    end
    combinedDcmdirS = newCombinedDcmdirS;
%     for i = 2:length(patNameC)
%         for j = 1:length(dcmdirS.(patNameC{i}).STUDY.SERIES)
%             combinedDcmdirS.STUDY.SERIES(end+1) = dcmdirS.(patNameC{i}).STUDY.SERIES(j);
%         end
%     end
for i = 2:length(patNameC)
    for  k = 1:length(dcmdirS.(patNameC{i}).STUDY)
        for j = 1:length(dcmdirS.(patNameC{i}).STUDY(k).SERIES)
            combinedDcmdirS.STUDY.SERIES(end+1) = dcmdirS.(patNameC{i}).STUDY(k).SERIES(j);
        end
    end
end
% Pass the java dicom structures to function to create CERR plan
    planC = dcmdir2planC(combinedDcmdirS);
else
    % Pass the java dicom structures to function to create CERR plan
    planC = dcmdir2planC(dcmdirS.(selected)); %wy
end
 
%-------------Store CERR version number---------------%
 
toc;
pause(0.05)
save_planC_QL(planC,'passed',bz2FileName);

