function [doseBinsV, volsHistV] = getDVHMatrix_QL(StuctureMask, Dose, numBins)
%"getDVHMatrix"
%   Wrapper function to get DVH data.  Use this to get the doseV/volsV
%   vectors used to represent DVH data for a structNum and doseNum.  If the
%   data has already been calculated and stored in the plan, the stored
%   data is returned, else new data is calculated using getDVH.
%
%   Matching the requested structNum and doseNum to existing DVHs is done
%   using the structure's name, and the dose index.  If a stored DVH has a
%   dose index that is no longer accurate, incorrect data may be returned.
%
%   planC is an output argument in order to save the calculated DVH data
%   for future use.
%
% Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team.
%
% This file is part of The Computational Environment for Radiotherapy Research (CERR).
%
% CERR development has been led by:  Aditya Apte, Divya Khullar, James Alaly, and Joseph O. Deasy.
%
% CERR has been financially supported by the US National Institutes of Health under multiple grants.
%
% CERR is distributed under the terms of the Lesser GNU Public License.
%
%     This version of CERR is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
% CERR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CERR.  If not, see <http://www.gnu.org/licenses/>.
%
%Usage:
%   [planC, doseBinsV, volsHistV] = getDVHMatrix(planC, structNum, doseNum)

%Find all DVHs for structure structNum on dose doseNum.
dosesV = Dose(StuctureMask==1);
volsV = 0.25^3*ones(numel(dosesV),1);
if(max(dosesV)>0)
    binWidth = max(dosesV)./numBins;
    %Histogram the volumes by dose.
    [doseBinsV, volsHistV] = doseHist(dosesV, volsV, binWidth);
else
    doseBinsV = 0*ones(numBins,1);
    volsHistV = volsV(1)*ones(numBins,1);
end
end