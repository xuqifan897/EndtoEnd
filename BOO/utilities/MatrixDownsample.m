function M = MatrixDownsample(Mhr,thresh)

[Ih,Jh,Vh] = find(Mhr);
[sizeI, sizeJ] = size(Mhr);
x = rand(size(Mhr,2),1);
Mhrx = Mhr*x;
whoMhr = whos('Mhr');
clearvars Mhr

cutoff = max(Vh)*thresh;
IndV = find(Vh>=cutoff);
Il = Ih(IndV);
Jl = Jh(IndV);
Vl = Vh(IndV);
clearvars Vh Ih Jh

M = sparse(Il,Jl,Vl,sizeI,sizeJ);
whoM = whos('M');
SizeScale = whoMhr.bytes/whoM.bytes;
fprintf(['\nMatrix downsampled by ' num2str(SizeScale) ' times\n']);

Mx = M*x;
fprintf(['\nError in M*x is ' num2str(norm(Mx-Mhrx)/norm(Mhrx)) ' \n']);

