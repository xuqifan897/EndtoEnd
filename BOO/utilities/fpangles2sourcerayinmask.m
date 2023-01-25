function [Pos, dirXray] = fpangles2sourcerayinmask(angles,isocenter,lengthinpixels)

phi = angles(:,1)/180*pi;
theta = angles(:,2)/180*pi;

x1 = sin(theta).*cos(phi);
y1 = sin(theta).*sin(phi);
z1 = cos(theta);

x = -y1;
y = -x1;
z = z1;

Pos = [x,y,z]*lengthinpixels+isocenter;
dirXray = [x,y,z];