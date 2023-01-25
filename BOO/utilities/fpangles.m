function theta=fpangles()
angres=6;
N=round(180/angres); deltaphi=180/N;m=0;
i=0; theta=0; alldose=0;
for thetay=0:deltaphi:180
    for thetax = 0:abs(round(360/(2*pi*cos((thetay-90)/180*pi)/(deltaphi/180*pi)))):359
        i=i+1;theta(i,1)=thetax;
        theta(i,2)=thetay;
    end
end

% for i=1:360
% theta(i,1:2)=[i-1,90];
% end