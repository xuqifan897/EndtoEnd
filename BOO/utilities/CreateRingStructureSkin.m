function [T_ringstruct, T_skin] = CreateRingStructureSkin(PTV, BODY, varargin)

radius=2*(3*sum(PTV(:))/4/pi)^(1/3);
if(radius>20)
    radius = 20;
end

if(isempty(varargin))
    gap = 3;
else
    gap = varargin{1};
end

[XS,YS,ZS]=meshgrid(linspace(-radius, radius, 2*radius),linspace(-radius, radius, 2*radius),linspace(-radius, radius, 2*radius));
binarysphere=(XS.^2 + YS.^2 + ZS.^2)<=radius^2;
T_ringstruct=logical(imdilate(PTV,binarysphere)-imdilate(PTV,ones(gap,gap,gap)));
T_ringstruct(BODY==0)=0;

T_skin=logical(BODY-imerode(BODY,ones(3,3,3)));
