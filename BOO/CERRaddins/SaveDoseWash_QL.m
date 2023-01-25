function [EntireImg,Imgs,ImgsInit,Masks] = SaveDoseWash_QL(patInfo,figuresFolder,figureName,DoseInd,patientName,FigureNum,ImageSize,varargin)
% Automatic generation of DoseWash figures including transverse, sagittal,
% and coronal directions
% Note that some parameters may need to be adjusted for proper use,
% including BWthresh,WidthRange,HeightRange1,HeightRange2

global planC stateS
mkdir([ figuresFolder '\temp\' ])
screensize = get( groot, 'Screensize' );
figure(FigureNum);set(gcf,'pos',screensize);
thresholdValue = 20;
SmoothFactor = 100;
margin = [10,10];
BWthresh = 50000;
WidthRange = [1/6,5/7];
HeightRange1 = [1/9,7/9];
HeightRange2 = [3.5/9,8/9];
pausetime = 0.5;

% load('D:\datatest\patient\patInfo.mat')
ii=find(strcmp({patInfo.Name},patientName));
StructureInfo = patInfo(ii).StructureInfo;
BODY = (imdilate((StructureInfo(2).Mask|StructureInfo(1).Mask),ones(3,3,3))*2000);
% BODY = (StructureInfo(2).Mask|StructureInfo(1).Mask)*2000;
addDoseToGui_dvo(BODY,'body')
planNum = size(planC{1,9},2);
sliceCallBack_QL('SELECTDOSE', planNum);

CERRStruct = patInfo(ii).CERRStruct;
CTwindow = [2000,2001];
doseDisplay = [1999,2000];
colorbarRange = [0,100];
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
sliceCallBack_QL('SELECTAXISVIEW',hAxis,'transverse',coordInd)

ImgDim = {'transverse','sagittal','coronal'};
Masks = cell(3,1);
figure(FigureNum);set(gcf,'pos',screensize);

for jj = 1:numel(ImgDim)
    if(~exist([figuresFolder '\temp\' figureName '_body_' ImgDim{jj} '.png'],'file'))
        sliceCallBack_QL('SELECTAXISVIEW',hAxis,ImgDim{jj},coordInd)
        pause(pausetime)
        robo = java.awt.Robot;
        t = java.awt.Toolkit.getDefaultToolkit();
        rectangle = java.awt.Rectangle(t.getScreenSize());
        image = robo.createScreenCapture(rectangle);
        filehandle = java.io.File([figuresFolder '\temp\' figureName '_body_' ImgDim{jj} '.png']);
        javax.imageio.ImageIO.write(image,'png',filehandle);
    end
    eval(['Img = rgb2gray(imread(''' figuresFolder '\temp\' figureName '_body_' ImgDim{jj} '.png''));']);
    Mask = Img>thresholdValue;
    [ImgH,ImgW] = size(Img);
    Ind = zeros(size(Mask));
    if jj == 1
        HeightRange = HeightRange1;
    else
        HeightRange = HeightRange2;
    end
    Ind(ImgH*HeightRange(1):ImgH*HeightRange(2),ImgW*WidthRange(1):ImgW*WidthRange(2)) = 1;
    Mask = Mask.*Ind;
    Mask = imerode(imdilate(Mask,ones([SmoothFactor,SmoothFactor])),ones([SmoothFactor,SmoothFactor]));
    Mask = imfill(Mask,'holes');
    Mask = bwareaopen(Mask,BWthresh);
    Masks{jj} = Mask ;
end

planC{1,9}(planNum) = [];

Imgs = cell(ImageSize(1),ImageSize(2)*3);

for mm = 1:length(DoseInd)
    ChangeCERRdoseWash_QL(patientName,patInfo,DoseInd(mm))
    for jj = 1:numel(ImgDim)
        if(~exist([figuresFolder '\temp\' figureName '_doseWash_' ImgDim{jj} '_doseInd_' num2str(mm) '.png'],'file'))
            sliceCallBack_QL('SELECTAXISVIEW',hAxis,ImgDim{jj},coordInd)
            pause(pausetime)
            robo = java.awt.Robot;
            t = java.awt.Toolkit.getDefaultToolkit();
            rectangle = java.awt.Rectangle(t.getScreenSize());
            image = robo.createScreenCapture(rectangle);
            filehandle = java.io.File([figuresFolder '\temp\' figureName '_doseWash_' ImgDim{jj} '_doseInd_' num2str(mm) '.png']);
            javax.imageio.ImageIO.write(image,'png',filehandle);
        end
        
        eval(['Img = imread(''' figuresFolder '\temp\' figureName '_doseWash_' ImgDim{jj} '_doseInd_' num2str(mm) '.png'');']);
        RowInd = ceil(mm/ImageSize(2));
        ColInd = mm - (RowInd-1)*ImageSize(2);
        Imgs{RowInd,(ColInd-1)*3+jj} = CropMaskImg_tight(Img, Masks{jj});
        ImgsInit{RowInd,(ColInd-1)*3+jj} = Img;
    end
end

EntireImg = PutImgTogether(Imgs,margin);
[nr,nc,~] = size(EntireImg);
ImageShowSize(1,2) = screensize(4);
ImageShowSize(1,1) = screensize(4)*nc/nr;
figure('pos',[1,1,ImageShowSize]);imagesc(EntireImg);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
axis off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
savefig(fig,[figuresFolder figureName '_doseWash.fig']);
print(fig, '-dtiff', [figuresFolder figureName '_doseWash.png'],'-r300');

end

