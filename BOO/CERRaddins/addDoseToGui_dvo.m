function addDoseToGui_dvo(bodyDoseArr,planName)

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
xVals = xoffset - (sizeDim2*grid2)/2 ;
yVals = yoffset + (sizeDim1*grid1)/2;
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