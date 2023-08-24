function Ratio = addDoseToGui_Normalize_Move_QL(bodyDoseArr,planName, presDose, PTVnum,xoff,yoff)
% Normalize the dose for the plan (# planNo) to presDose
% Move the dose when there is a misalignment
% The xoff and yoff is tuned by hand with the following steps:
% 1. set xoff and yoff: e.g. xoff = 0; yoff = 1;
% 2. add PTV to CERR: addDoseToGui_Move_QL(PTV,planName,xoff,yoff)
% 3. see if the PTV match with the PTV shown in CERR
% where PTV is the PTV mask used in the optimization

addDoseToGui_Move_QL(bodyDoseArr,planName,xoff,yoff)

global planC

planNum = size(planC{1,9},2);

doseNum = planNum;

[dosesV, volsV]=getDVH(PTVnum, doseNum, planC);
[doseBinsV, volsHistV] = doseHist(dosesV, volsV, 0.2);
cumVolsV = cumsum(volsHistV);
cumVols2V  =(cumVolsV(end)-cumVolsV)/ cumVolsV(end);  %cumVolsV is the cumulative volume lt that corresponding dose
Ind1 = find(cumVols2V<0.95, 1 );
Ind2 = Ind1 - 1;

D95 = (doseBinsV(Ind1)-doseBinsV(Ind2))/(cumVols2V(Ind1)-cumVols2V(Ind2))*(0.95-cumVols2V(Ind1))+doseBinsV(Ind1);

Ratio = presDose/D95;
dose_norm = bodyDoseArr*Ratio;

planC{1,9}(planNum) = [];

addDoseToGui_Move_QL(dose_norm,planName,xoff,yoff)

end

