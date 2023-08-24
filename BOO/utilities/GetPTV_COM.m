function [CenterOfMass] = GetPTV_COM(PTV)

ImageX=sum(sum(PTV,3),2);
MassSum=sum(ImageX);
Xc=0;Yc=0;Zc=0;
for i = 1:length(ImageX)
    Xc = Xc + i*ImageX(i)/MassSum;
end
Xc = floor(Xc+0.5);

ImageY=sum(sum(PTV,3),1);
for i = 1:length(ImageY)
    Yc = Yc + i*ImageY(i)/MassSum;
end
Yc = floor(Yc+0.5);

ImageZ=sum(sum(PTV,2),1);
for i = 1:length(ImageZ)
    Zc = Zc + i*ImageZ(i)/MassSum;
end
Zc = floor(Zc+0.5);

CenterOfMass = [Xc, Yc, Zc];


end
