function out = draw_beammask_QL(fpang,BODY,PTV)
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

bin=500;
max_r = 65;
r=0:max_r/(bin-1):max_r;

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

Az = -88;
El = 34;
figure('position',[-74 415 1261 844]);
hold on

ds=BODY; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','red','EdgeColor','none');
isonormals(ds,pct);
view(Az,El);camlight left; camlight right; lighting gouraud;

ds=PTV; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','green','EdgeColor','none');
isonormals(ds,pct); view(Az,El);camlight left; camlight right; lighting gouraud;

ds=A; d50=max(max(max(ds)))*0.5;pct=patch(isosurface(ds,d50/100),'FaceColor','blue','EdgeColor','none');
isonormals(ds,pct); view(Az,El);camlight left; camlight right; lighting gouraud;

axis([-20 201 -20 180 0 173]);
axis equal
axis tight
axis off
alpha(0.5)

set(gcf,'color','w');


clear ds d50 pct;
