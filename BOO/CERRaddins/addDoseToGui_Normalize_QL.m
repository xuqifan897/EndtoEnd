function scale = addDoseToGui_Normalize_QL(bodyDoseArr,planName, presDose, PTVnum)
% Normalize the dose for the plan (# planNo) to presDose

addDoseToGui_dvo(bodyDoseArr,planName);
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

scale = presDose/D95;
dose_norm = bodyDoseArr*presDose/D95;

planC{1,9}(planNum) = [];

addDoseToGui_dvo(dose_norm,planName);

end

