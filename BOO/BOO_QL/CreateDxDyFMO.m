function [Dx,Dy] = CreateDxDyFMO (BeamletLog0)

[Nx,Ny,Nbeams] = size(BeamletLog0);

Dxs=sparse(-eye(Nx));
Dxs_1=sparse(diag(ones(Nx-1,1),1));
Dxs(Nx,Nx)=0;
Dxs=Dxs+Dxs_1;
IN=sparse(eye(Ny));
Dx0=kron(IN,Dxs);

Dys=sparse(-eye(Ny));
Dys_1=ones(Ny-1,1);
Dys_1=sparse(diag(Dys_1,1));
Dys(Ny,Ny)=0;
Dys=Dys+Dys_1;
IM=sparse(eye(Nx));
Dy0=kron(Dys,IM);%x

Dx=kron(eye(Nbeams),Dx0);
Dy=kron(eye(Nbeams),Dy0);

Dx(:,BeamletLog0==0)=[];
Dy(:,BeamletLog0==0)=[];

Dxsum = sum(Dx,2);
Dysum = sum(Dy,2);

Dx(Dxsum~=0,:) = [];
Dy(Dysum~=0,:) = [];

Dxsum = sum(abs(Dx),2);
Dysum = sum(abs(Dy),2);

Dx(Dxsum==0,:) = [];
Dy(Dysum==0,:) = [];
