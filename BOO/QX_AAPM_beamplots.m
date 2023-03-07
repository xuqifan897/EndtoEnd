% This is to draw the plots of the candidate beams for AAPM abstract submission

addpath(genpath('BOO_QL'), '-end');
addpath(genpath('CERR2016'), '-end');
addpath(genpath('CERRaddins'), '-end');
addpath(genpath('utilities'), '-end');
addpath(genpath('beamlogs'), '-end');

globalFolder = '/data/datasets/UCLAPatients';
expFolder = fullfile(globalFolder, 'experiment');
visFolder = fullfile(globalFolder, 'visNew');
numPatients = 8;

for i = 1:numPatients
    patientName = ['patient', num2str(i)];
    expPatFolder = fullfile(expFolder, patientName);
    optFolder = fullfile(expPatFolder, 'optimizePTVcropped');

    % load result
    pattern = fullfile(optFolder, 'BOOresult*.mat');
    files = dir(pattern);
    trailNOs = [];
    for j = 1:length(files)
        digitFlag = isstrprop(files(j).name, 'digit');
        trailNO = str2num(files(j).name(digitFlag));
        trailNOs(end+1) = trailNO;
    end
    trailNO = max(trailNOs);
    BOOResultFile = fullfile(optFolder, ['BOOresult', num2str(trailNO), '.mat']);
    
    load(BOOResultFile, 'BOOresult');
    beamfpangles = BOOresult.params.beamfpangles;
    activeBeams = BOOresult.activeBeams;
    PTV = BOOresult.StructureInfo(1).Mask;
    BODY = BOOresult.StructureInfo(2).Mask;
     draw_beammask_QL(beamfpangles(activeBeams, :), BODY, PTV);

     outFile = fullfile(visFolder, patientName, 'BOObeams.png');
     saveas(gcf, outFile);
     close(gcf);
end