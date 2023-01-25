%load segment PTV BODY
load('D:\dicom\4pi\SPNAs1\segment.mat','BODY','PTV');
load('D:\fmo4pi\optimize\SPNAs1_50gy\beamletIntOutput20.mat');
ang=fpangles();
selangles=ang(beamnum,:);
fpang=selangles;
%angles=fpangles();fpang=angles(beamnum,:);
theta=fpang(:,1)*pi/180;
phi=fpang(:,2)*pi/180;

x=sum(squeeze(sum(PTV,3)),1);
locationcount=1:length(x);
Iy=round(sum(locationcount.*x)/sum(x));

x=sum(squeeze(sum(PTV,3)),2);
locationcount=1:length(x);
locationcount=locationcount';
Ix=round(sum(locationcount.*x)/sum(x));

x=sum(squeeze(sum(PTV,1)),1);
locationcount=1:length(x);
Iz=round(sum(locationcount.*x)/sum(x));
[sx,sy,sz]=size(PTV);

[X,Y,Z] = meshgrid(1-Iy:sy-Iy,1-Ix:sx-Ix,1-Iz:sz-Iz);
xindex=Y(:,1,1);
yindex=X(1,:,1);
zindex=squeeze(Z(1,1,:));

bin=round(norm(size(PTV)));
r=0:bin/(bin-1):bin;

A=zeros(size(PTV));
for p=1:length(theta)
    [x,y,z]=sph2cart(ones(1,bin)*(3*pi/2-theta(p)),ones(1,bin)*(pi/2-phi(p)),r);
    
    for i=1:bin
        if x(i)> max(xindex) || x(i) < min(xindex) || y(i)> max(yindex) || y(i) < min(yindex) || z(i)> max(zindex) || z(i) < min(zindex)
            continue
        end
        [val locx]=min(abs(x(i)-xindex));
        [val locy]=min(abs(y(i)-yindex));
        [val locz]=min(abs(z(i)-zindex));
        A(locx,locy,locz)=1;
    end
end

figure();
hold on
ds=BODY; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','red','EdgeColor','none');view(92,0);camlight left; camlight; lighting gouraud;
ds=PTV; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','green','EdgeColor','none');view(92,0);camlight left; camlight; lighting gouraud;
ds=A; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','blue','EdgeColor','none');view(92,0);camlight left; camlight; lighting gouraud;
%ds=CenterClarge; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','yellow','EdgeColor','none');view(92,0);camlight left; camlight; lighting gouraud;
axis([0 size(PTV,1) 0 size(PTV,2) 0 size(PTV,3)]);
axis equal
alpha(0.5)

clear ds d50 pct;
