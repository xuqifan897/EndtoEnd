function dose0new = getDoseFromCERR(DoseInd, RefDoseInd)
% export dose from CERR, based on the dimension of a reference dose

global planC

x0 = planC{1, 9}(DoseInd).horizontalGridInterval;
y0 = planC{1, 9}(DoseInd).verticalGridInterval;
z0 = mean(diff(planC{1, 9}(DoseInd).zValues));

dosesize0 = size(planC{1, 9}(DoseInd).doseArray);
dose0 = planC{1, 9}(DoseInd).doseArray;

x1 = planC{1, 9}(RefDoseInd).horizontalGridInterval;
y1 = planC{1, 9}(RefDoseInd).verticalGridInterval;
z1 = mean(diff(planC{1, 9}(RefDoseInd).zValues));

dosesize1 = size(planC{1, 9}(RefDoseInd).doseArray);
dose1 = planC{1, 9}(RefDoseInd).doseArray;

newsize = round(dosesize0./[x1/x0,y1/y0,z1/z0]);
newsize = min(newsize,dosesize1);
dose0 = imresize3(dose0,newsize);

xinit0 = planC{1, 9}(DoseInd).coord1OFFirstPoint;
yinit0 = planC{1, 9}(DoseInd).coord2OFFirstPoint;
zinit0 = planC{1, 9}(DoseInd).zValues(1);
xinit1 = planC{1, 9}(RefDoseInd).coord1OFFirstPoint;
yinit1 = planC{1, 9}(RefDoseInd).coord2OFFirstPoint;
zinit1 = planC{1, 9}(RefDoseInd).zValues(1);

sx = 1-floor((xinit1-xinit0)/x1);
sy = 1+floor((yinit1-yinit0)/x1);
sz = 1-floor((zinit1-zinit0)/z1);

dose0new = zeros(dosesize1);
dose0new(sy:sy+newsize(1)-1,sx:sx+newsize(2)-1,sz:sz+newsize(3)-1)=dose0;
figure;imshow3D([dose1 dose0new dose1-dose0new])

dose0new = flip(dose0new,3);

