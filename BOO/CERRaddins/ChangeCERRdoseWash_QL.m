function ChangeCERRdoseWash_QL(patientName,patInfo,varargin)
% This function set up DoseWash in CERR for publication
% Parameters in patInfo may need manual set up
% An example of patInfo can be found in 'patInfo.mat'
% Refer to 'SaveDoseWash_QL.m' for automatic generation of DoseWash figures

% load('patInfo.mat')
ii=find(strcmp({patInfo.Name},patientName));

CERRStruct = patInfo(ii).CERRStruct;
CTwindow = patInfo(ii).CTwindow;
doseDisplay = patInfo(ii).doseDisplay;
colorbarRange = patInfo(ii).colorbarRange;
StructureInfo = patInfo(ii).StructureInfo;

global planC stateS

if(~isempty(varargin))
    planNum = varargin{1};
else
    planNum = size(planC{1,9},2);
end

sliceCallBack_QL('SELECTDOSE', planNum);
sliceCallBack_QL('DOSETOGGLE','on')
sliceCallBack_QL('TOGGLESINGLESTRUCT',CERRStruct)
sliceCallBack_QL('PLANELOCATORTOGGLE','off')
sliceCallBack_QL('CTWINDOW',CTwindow) % Dose/CT
sliceCallBack_QL('SLIDERTRANSALPHA',0.8) % Dose/CT
stateS.optS.staticColorbar = 1;
stateS.doseDisplayRange = doseDisplay;
stateS.colorbarRange = colorbarRange;
stateS.doseDisplayChanged   = 1;
stateS.colorbarChanged     = 1;
CERRRefresh
hAxis = stateS.handle.CERRAxis(1);
if(isempty(patInfo(ii).coordInd))
    coordInd = getPTVcoordInd(StructureInfo(1).Mask);
else
    coordInd = patInfo(ii).coordInd;
end
sliceCallBack_QL('SELECTAXISVIEW',hAxis,'transverse',coordInd) %% 3: 'transverse' 1: 'sagittal' 2: 'coronal'
end