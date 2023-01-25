%% This is an example of how to use the CERR functions in the CERRaddins folder

%% Useful functions
% Import all dicom files in DicomPath and its subfolders
% save planC in OutputFileName
CERRImportDicom_QL(DicomPath,bz2FileName); % bz2FileName = 'patientName.mat.bz2.mat';

% Open CERR viewer
CERR('CERRSLICEVIEWER')

% Load planC
sliceCallBack_QL('OPENNEWPLANC', bz2FileName);  % bz2FileName = 'patientName.mat.bz2.mat';

% Select dose
sliceCallBack_QL('SELECTDOSE', num2str(doseInd));

% View dose toggle
sliceCallBack_QL('DOSETOGGLE','on')

% Plane Locator toggle
sliceCallBack_QL('PLANELOCATORTOGGLE','off')

% Change colorbar and dose display
stateS.doseDisplayRange = [4,43];
stateS.colorbarRange = [0,42];
CERRRefresh

%% Not clear
CERRAxisMenu('SET_VIEW')
CERRAxisMenu('UPDATE_MENU')
hAxis = get(gcbo, 'userdata');
sliceCallBack_QL('SELECTAXISVIEW',hAxis,'coronal')
sliceCallBack_QL('AXISCLICKED')

