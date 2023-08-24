function [scale,D95] = plotDVH_QL(DoseInfo, strNum, StructureInfo, numBins, varargin)
% plot DVH for selected structures in StructureInfo of all plans in DoseInfo
% strNum: the index of selected structures in StructureInfo
% must have CERR opened
% the colors were specified in DVHcolormap.mat

if(length(varargin)>=1)
    Normalize = varargin{1};
else
    Normalize = 1;
end
global planC
global stateS

cum_diff_string = 'CUMU';

Dose = {DoseInfo.Data};
presDose = StructureInfo(strNum(1)).IdealDose;
numPlan = numel(Dose);
numStruct = numel(strNum);

indexS = planC{end};
% clear previous DVH cache
planC{indexS.DVH}(1:end) = [];

% load DVH information into planC(indexS.DVH)
for iPlan=1:numPlan
    for iStruct=1:numStruct
        DVHIndex = length(planC{indexS.DVH}) + 1;
        [doseBinsV, volsHistV] = getDVHMatrix_QL(StructureInfo(strNum(iStruct)).Mask==1, Dose{iPlan}, numBins);
        if(iStruct==1)
            cumVolsV = cumsum(volsHistV);
            cumVols2V  =(cumVolsV(end)-cumVolsV)/ cumVolsV(end);  %cumVolsV is the cumulative volume lt that corresponding dose
            Ind1 = find(cumVols2V<0.95, 1 );
            Ind2 = Ind1 - 1;
            D95 = (doseBinsV(Ind1)-doseBinsV(Ind2))/(cumVols2V(Ind1)-cumVols2V(Ind2))*(0.95-cumVols2V(Ind1))+doseBinsV(Ind1);
            if(Normalize)
                scale(iPlan) = presDose/D95;
            else
                scale(iPlan) = 1;
            end
            Dose{iPlan} = Dose{iPlan}*scale(iPlan);
            
            [doseBinsV, volsHistV] = getDVHMatrix_QL(StructureInfo(strNum(iStruct)).Mask, Dose{iPlan}, numBins);
        end
        planC = saveDVHMatrix(DVHIndex, doseBinsV, volsHistV, planC);
        %Create a new DVH element.
        structureName = StructureInfo(strNum(iStruct)).Name;
        planC{indexS.DVH}(DVHIndex).structureName = structureName;
        planC{indexS.DVH}(DVHIndex).volumeType      = 'absolute';
        planC{indexS.DVH}(DVHIndex).doseType        = 'absolute';
        planC{indexS.DVH}(DVHIndex).doseIndex       = iPlan;
        planC{indexS.DVH}(DVHIndex).assocStrUID=strNum(iStruct);
    end
end

legend_str = '';
legend_string = '';
legend_string_abs = '';
gridFlag = 0;

%Prepare a figure.
h = figure('tag', 'DVHPlot', 'doublebuffer', 'on');
pos = get(h,'position');
nDVHFigs = length(findobj('tag', 'DVHPlot'));
set(h,'position',[pos(1)*(1 - 0.05*nDVHFigs),pos(2)*(1 - 0.05*nDVHFigs),pos(3),pos(4)])
set(gcf,'units','normalized','outerposition',[0 0 1 1])

try
    stateS.handle.DVHPlot = [stateS.handle.DVHPlot h];
catch
    stateS.handle.DVHPlot = h;
end

set(h, 'userdata', []);
uimenu(h, 'label', 'Expand Options', 'callback',['plotDVHCallback(''EXPANDEDVIEW'')'],'interruptible','on');
set(h,'name',['DVH plot: ' stateS.CERRFile])
color2V = get(stateS.handle.CERRSliceViewer, 'Color');
set(h,'color', color2V);
set(h,'numbertitle','off')
hAxis = axes('parent', h);

ylabel('Fractional volume', 'parent', hAxis)
title('Dose volume histograms', 'parent', hAxis)

grid(hAxis, 'off');

hOld = findobj('tag','CERRAbsDVHPlot');
gridSetting = 'on'; %default
try
    figure(hOld(1))
    gridSetting = get(gca,'xgrid');
    delete(hOld)
end

volV = logical(ones(numPlan*numStruct,1))';
surfV   = logical(zeros(numPlan*numStruct,1))';
avgV    = logical(zeros(numPlan*numStruct,1))';
absV    = logical(zeros(numPlan*numStruct,1))';

%Iterate over the volV, surfV DVH lists, calculate if flagged.
skipDVH = 0;

if length(volV)==1
    dispLegend = 0;
else
    dispLegend = 1;
end

oldUnits = '';

count = 1;
hold on
for i = 1 : length(volV)
    
    DVHNum  = i;
    doseSet = [];
    str = planC{indexS.DVH}(i).structureName;
    structNum = planC{indexS.DVH}(i).assocStrUID;
    
    if structNum ~= 0
%         if isempty(planC{indexS.structures}(structNum).rasterSegments)
%             warning(['No scan segments stored for structure ' num2str(structNum) ])
%         end
    end
    sNames = {planC{indexS.DVH}.structureName};
    sameName = find(strcmpi(sNames, str) & volV);
    flagLSS = find(sameName == i);
    if isempty(doseSet) || isempty(structNum)
        units = planC{indexS.DVH}(i).doseUnits;
    else
        units = getDoseUnitsStr(doseSet,planC);
    end
    if isempty(oldUnits)
        xlabel(['Dose ' units], 'parent', hAxis)
        oldUnits = units;
    elseif ~isempty(oldUnits) && ~strcmpi(units,oldUnits)
        errordlg('Dose units must be same for all DVHs.','Dose Units','modal')
        return;
    end
    
    %Draw the volume DVH.
    %[planC] = showDVH(hAxis, DVHNum, doseSet, gridSetting, flagLSS, absV(i), gridFlag);
    doseStat = showDVH(hAxis, DVHNum, doseSet, gridSetting, flagLSS, absV(i), gridFlag, cum_diff_string);
    if isempty(doseStat)
        return
    end
    drawnow
    name = planC{indexS.DVH}(DVHNum).structureName;
    if ~isempty(planC{indexS.DVH}(DVHNum).fractionIDOfOrigin)
        name = [name, ' (', num2str(planC{indexS.DVH}(DVHNum).fractionIDOfOrigin), ')'];
    end
    name = regexprep(name,'_','-');
    if dispLegend % if only one DVH is displayed do not show legend
        if isempty(legend_string)
            legend_string = name;
        else
            legend_string = char(legend_string,name); %strvcat(legend_string,name);
        end
        
        if absV(i) && isempty(legend_string_abs)
            legend_string_abs = name;
        elseif absV(i)
            legend_string_abs = char(legend_string_abs,name);
        end
        legend(hAxis, legend_string,'Location', 'NorthEastOutside');%adds legend to DVH plot
        absDVHAxis = findobj('tag', 'AbsDVHAxis');
        if ~isempty(absDVHAxis)
            legend(absDVHAxis,'Location', 'NorthEastOutside', legend_string_abs);%adds legend to DVH plot
        end
    end
end

planC{indexS.DVH}(1:end) = [];

A=get(legend(gca),'String');
for i=1:numStruct
    B{i}=A{i};
end

h=legend(B,'FontSize',30);
set(h,'Location','northeastoutside');
set(gca,'color','None');

hline = findobj(gcf, 'type', 'line');
load('DVHcolormap.mat');
% load('colormap.mat');

hline = flip(hline);
for i =1:numPlan
    for j = 1:numStruct
        %         hline((i-1)*numStruct+j).Color
        hline((i-1)*numStruct+j).Color = DVHcolormap(mod(j,length(DVHcolormap)-1)+1,:);
    end
end
set(hline,'LineWidth',4)
xlabel('Dose (Gy)','FontSize',30);
ylabel('Fractional volume','FontSize',30);
set(gca,'fontsize',30);

DoseName = [];
for i =1:numPlan
    tmp = strrep(DoseInfo(i).Name,'_',' ');
    if i==1
        tmp = ['Solid: ' tmp];
    elseif i==2
        tmp = ['Dotted: ' tmp];
    elseif i==3
        tmp = ['Dash dotted: ' tmp];
    elseif i==4
        tmp = ['Dashed: ' tmp];
    end
    DoseName = [DoseName tmp '; '];
end
title(DoseName,'FontSize',20);
grid off

end
