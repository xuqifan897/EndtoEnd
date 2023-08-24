function SaveStructureFileOnly(ROINameList, jsonFileName, PresDose)
% Save jason file containing structures for dose calculation 
NumStruct = numel(ROINameList);

for ii = 1:NumStruct
    disp(['Structure ' num2str(ii) ': ' ROINameList{ii}]);
end

PTVind = input('Input PTV index: ');
BODYind = input('Input BODY index: ');

OARname = {};
count = 1;
for ii = 1:NumStruct
    ROIName = ROINameList{ii};
    if(ii==PTVind)
        PTVname = ROIName;
    elseif(ii==BODYind)
        BODYname = ROIName;
    else
        OARname{count} = ROIName;
        count = count + 1;
    end
end

data = [];
data.prescription = PresDose;
data.ptv = PTVname;
data.oar = string([BODYname OARname])';
saveJSONfile(data, jsonFileName);

end


function saveJSONfile(data, jsonFileName)
% saves the values in the structure 'data' to a file in JSON format.
%
% Example:
%     data.name = 'chair';
%     data.color = 'pink';
%     data.metrics.height = 0.3;
%     data.metrics.width = 1.3;
%     saveJSONfile(data, 'out.json');
%
% Output 'out.json':
% {
% 	"name" : "chair",
% 	"color" : "pink",
% 	"metrics" : {
% 		"height" : 0.3,
% 		"width" : 1.3
% 		}
% 	}
%
    fid = fopen(jsonFileName,'w');
    %fid=1;
    writeElement(fid, data,'');
    fprintf(fid,'\n');
    fclose(fid);
end

% function buildJSON(fid,data)
%     namesOfFields = fieldnames(data);
%     numFields = length(namesOfFields);
%     
%     if isstruct(
%     fprintf(fid,'{ "%s" : [\n ',rootElementName);
%     for m = 1:length(data) - 1
%        fprintf(fid,'{  \n ');
%        writeSingleGene(fid,numFields, namesOfFields,dataname,data,m);
%        fprintf(fid,'},\n');
%     end
%     m= m+1;
%     fprintf(fid,'{  \n ');
%     writeSingleGene(fid, numFields, namesOfFields,dataname,data,m);
%     fprintf(fid,' }\n]\n');
%     
%    
%     end

function writeElement(fid, data,tabs)
    namesOfFields = fieldnames(data);
    numFields = length(namesOfFields);
    tabs = sprintf('%s\t',tabs);
    fprintf(fid,'{\n%s',tabs);
   

    for i = 1:numFields - 1
        currentField = namesOfFields{i};
        currentElementValue = eval(sprintf('data.%s',currentField));
        writeSingleElement(fid, currentField,currentElementValue,tabs);
        fprintf(fid,',\n%s',tabs);
    end
    if isempty(i)
        i=1;
    else
      i=i+1;
    end
      
    
    currentField = namesOfFields{i};
    currentElementValue = eval(sprintf('data.%s',currentField));
    writeSingleElement(fid, currentField,currentElementValue,tabs);
    fprintf(fid,'\n%s}',tabs);
end

function writeSingleElement(fid, currentField,currentElementValue,tabs)
    
        % if this is an array and not a string then iterate on every
        % element, if this is a single element write it
        if length(currentElementValue) > 1 && ~ischar(currentElementValue)
            fprintf(fid,' "%s" : [',currentField);
            for m = 1:length(currentElementValue)-1
                fprintf(fid,' "%s",',currentElementValue(m));
            end
            if isempty(m)
                m=1;
            else
              m=m+1;
            end
            fprintf(fid,' "%s"',currentElementValue(m));
            fprintf(fid,'%s]%s',tabs,tabs);
        elseif isstruct(currentElementValue)
            fprintf(fid,'"%s" : ',currentField);
            writeElement(fid, currentElementValue,tabs);
        elseif isnumeric(currentElementValue)
            fprintf(fid,'"%s" : %g' , currentField,currentElementValue);
        elseif isempty(currentElementValue)
            fprintf(fid,'"%s" : null' , currentField,currentElementValue);
        else %ischar or something else ...
            fprintf(fid,'"%s" : "%s"' , currentField,currentElementValue);
        end
end

% function writeSingleGene(fid,numFields, namesOfFields,dataname,data,m)
% 
%     for i = 1:numFields - 1
%         field = namesOfFields{i};
%         stmp = sprintf('%s(m).%s',dataname,field);
%         fprintf(fid,'"%s" : "%s",\n' , field,eval(stmp));
%     end
%     i=i+1;
%     field = namesOfFields{i};
%     stmp = sprintf('%s(m).%s',dataname,field);
%     fprintf(fid,'"%s" : "%s"\n',field,eval(stmp));
%         
%         
% end
