function [zendpos,zendmaskPos,zstartpos, zstartmaskPos] = GetBODYend_dimenless(BODY)

siz = size(BODY);
ImageZ=sum(sum(BODY,2),1);

zendpos = find(ImageZ,1,'first');
[I,J] = ind2sub(siz(1:2),find(BODY(:,:,zendpos)));
zendmaskPos = ([I,J,zendpos*ones(size(I))])';

zstartpos = find(ImageZ,1,'last');
[I,J] = ind2sub(siz(1:2),find(BODY(:,:,zstartpos)));
zstartmaskPos = ([I,J,zstartpos*ones(size(I))])';

