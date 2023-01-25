function addDoseToGui_Move_QL(bodyDoseArr,planName,xoff,yoff)
% Move the dose when there is a misalignment
% The xoff and yoff is tuned by hand with the following steps:
% 1. set xoff and yoff: e.g. xoff = 0; yoff = 1;
% 2. add PTV to CERR: addDoseToGui_Move_QL(PTV,planName,xoff,yoff)
% 3. see if the PTV match with the PTV shown in CERR
% where PTV is the PTV mask used in the optimization

global planC

fname = planName;
fullmaskValue = bodyDoseArr;

ds=fullmaskValue;ds0=ds;
sizex=size(planC{1,3}(1,1).scanArray,2);
sizey=size(planC{1,3}(1,1).scanArray,1);
sizet=size(planC{1,3}(1,1).scanArray,3);
xoffset=planC{1,3}(1,1).scanInfo(1,1).xOffset;
yoffset=planC{1,3}(1,1).scanInfo(1,1).yOffset;

sizeDim1=sizey-1;
sizeDim2=sizex-1;
sizeDim3=sizet-1;
grid1=planC{1,3}(1,1).scanInfo(1,1).grid1Units;
grid2=planC{1,3}(1,1).scanInfo(1,1).grid2Units;
grid3=planC{1,3}(1,1).scanInfo(1,1).voxelThickness;
xVals = xoffset - (sizeDim2*grid2)/2 + xoff;
yVals = yoffset + (sizeDim1*grid1)/2 + yoff;
ctunit1=grid1*sizeDim1/(size(fullmaskValue,1)-1);

regParamsS.horizontalGridInterval=ctunit1;
regParamsS.verticalGridInterval=-ctunit1;
regParamsS.coord1OFFirstPoint=xVals;
regParamsS.coord2OFFirstPoint=yVals;
ctunit3=grid3*sizeDim3/(size(fullmaskValue,3)-1);
regParamsS.zValues=planC{1,3}(1,1).scanInfo(1,1).zValue:ctunit3:(planC{1,3}(1,1).scanInfo(1,1).zValue+ctunit3*(size(fullmaskValue,3)-1));
assocScanNum = 1;indexS = planC{end};
assocScanUID = planC{indexS.scan}(assocScanNum).scanUID;
%    name=[beamconf,num];
 name=[fname];
%     name=[beamconf,num,'strat05'];
ds=flipdim(ds,3);
dose2CERR(ds,[],name,name,name,'' ,regParamsS, 'no', assocScanUID);

planNum = size(planC{1,9},2);
sliceCallBack_QL('SELECTDOSE', planNum);

end
