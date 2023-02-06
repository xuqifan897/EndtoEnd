% this code runs on CPU server, which requires large storage

addpath(genpath('BOO_QL'), '-end');
addpath(genpath('CERR2016'), '-end');
addpath(genpath('CERRaddins'), '-end');
addpath(genpath('utilities'), '-end');
addpath(genpath('beamlogs'), '-end');

globalFolder = '/data/datasets/UCLAPatients';
dataFolder = fullfile(globalFolder, 'anonymousDataNew');
expFolder = fullfile(globalFolder, 'experiment');
numPatients = 1;
thresh = 1e-06;
PrescriptionDose = 20;

% for idx = 1:numPatients
    patientName = ['patient', num2str(idx)];
    patExpFolder = fullfile(expFolder, patientName);
    optFolder = fullfile(patExpFolder, 'optimize');
    if ~ isfolder(optFolder)
        mkdir(optFolder)
    end
    h5file = fullfile(patExpFolder, 'Dose_Coefficients.h5');
    maskfile = fullfile(patExpFolder, 'Dose_Coefficients.mask');
    [M, dose_data, masks] = BuildDoseMatrix(h5file, maskfile, thresh);
    MorgFile = fullfile(patExpFolder, [patientName, '_M_original.mat']);
    save(MorgFile, 'M', 'dose_data', 'masks', '-v7.3');

    % save parameter files
    [StructureInfo, params] = InitIMRTparams(M, dose_data, masks, PrescriptionDose);
    ParamsNum = 0;
    paramsPath = fullfile(optFolder, ['params', num2str(ParamsNum), '_original.mat']);
    save(paramsPath, 'params');
    InfoNum = 0;
    infoPath = fullfile(optFolder, ['StructureInfo', num2str(InfoNum), '.mat']);
    save(infoPath, 'StructureInfo')
    %% Remove beams going through cutted off CT images
    PTV = StructureInfo(1).Mask;
    BODY = StructureInfo(2).Mask;
    figure;imshow3D([PTV, BODY],[])
    
    Isos = GetPTV_COM(PTV);
    [zendpos,zendmaskPos,zstartpos, zstartmaskPos] = GetBODYend_dimenless(BODY);
    
    beamfpangles = params.beamfpangles;
    [srcpos, dirXray] = fpangles2sourcerayinmask(beamfpangles,Isos,3000);
    
    dr2 = 50;
    numbeams = size(beamfpangles,1);
    badbeams = zeros(numbeams,1);
    for ii = 1:size(beamfpangles,1)
        isrcpos = srcpos(ii,:);
        idirXray = dirXray(ii,:);
        if(sign((zendpos-Isos(3))*(isrcpos(3) - Isos(3)))==1)
            if((zendpos-Isos(3))/idirXray(3)<0)
                error('Wrong direction!')
            end
            endpos = Isos + (zendpos-Isos(3))/idirXray(3)*idirXray;
            if(~isempty(zendmaskPos))
                if(find(sum((endpos'-zendmaskPos).^2,1)<dr2))
                    badbeams(ii) = 1;
                end
            end
        elseif(sign((zstartpos-Isos(3))*(isrcpos(3) - Isos(3)))==1)
            if((zstartpos-Isos(3))/idirXray(3)<0)
                error('Wrong direction!')
            end
            endpos = Isos + (zstartpos-Isos(3))/idirXray(3)*idirXray;
            if(~isempty(zstartmaskPos))
                if(any(sum((endpos'-zstartmaskPos).^2,1)<dr2))
                    badbeams(ii) = 1;
                end
            end
        end
    end
    
    % draw_beammask_QL(params.beamfpangles(1,:),BODY,PTV);
    draw_beammask_QL(params.beamfpangles(badbeams==1,:),BODY,PTV);
    draw_beammask_QL(params.beamfpangles(badbeams==0,:),BODY,PTV);
    
    goodbeamind = find(badbeams==0);
    newparams = params;
    newparams.beamWeightsInit = params.beamWeightsInit(goodbeamind);
    newparams.beamSizes = params.beamSizes(goodbeamind);
    newparams.BeamletLog0 = params.BeamletLog0(:,:,goodbeamind);
    newparams.beamVarianIEC = params.beamVarianIEC(goodbeamind,:);
    newparams.beamfpangles = params.beamfpangles(goodbeamind,:);
    
    BeamletLog0 = params.BeamletLog0;
    BeamletLog0Ind = BeamletLog0;
    BeamletLog0Ind(BeamletLog0==1) = 1:nnz(BeamletLog0);
    newbeamletind = BeamletLog0Ind(:,:,goodbeamind);
    newbeamletind = newbeamletind(newbeamletind>0);
    
    if(sum(newparams.beamSizes)~=numel(newbeamletind))
        error('Dimension error in removing beams!')
    end
    
    M = M(:,newbeamletind);
    save(fullfile(patExpFolder,[patientName '_M.mat']),'M','-v7.3');
    
    params = newparams;
    save(fullfile(optFolder,['params' num2str(ParamsNum) '.mat']),'params');
% end
exit;