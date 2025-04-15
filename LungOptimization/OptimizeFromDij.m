%% Read csv file
datafile = readtable('SingleLesionAllEQD2.csv');
patients = readlines('test.txt');
doseParentDir = "E:\matRadData\SingleLesionLungDosesPreprocessed\";
doseSaveDir = "E:\matRadData\SingleLesionPhysicalDosesOptimizedFromPreprocessed\";

%%
for i = 1:numel(patients)
    close all
    patient = patients(i);
    %runInversePlanning(patient, datafile, fullfile(doseParentDir, "GAN" + patient), fullfile(doseSaveDir, "GAN" + patient + ".mat"));
    %runInversePlanning(patient, datafile, fullfile(doseParentDir, "Unet" + patient),  fullfile(doseSaveDir, "Unet" + patient + ".mat"));
    runInversePlanning(patient, datafile, fullfile(doseParentDir, "HDUnet" + patient), fullfile(doseSaveDir, "HDUnet" + patient + ".mat"))
    %runInversePlanning(patient, datafile,  fullfile(doseParentDir, patient), fullfile(doseSaveDir, patient + ".mat"))
end

%% 
%function resultGUI3 = runInversePlanning(patient, datafile, pathToPredDose, result_save_path, dose_save_path)
function runInversePlanning(patient, datafile, pathToPredDose, dose_save_path)
    patientParentDir = "E:\matRadData\SingleLesionLungPatients\";
    dijDir = "E:\\matRadData\dijPBK3";
    
    tic;
    matRad_rc
    load(fullfile(patientParentDir, patient));
    clear resultGUI %Don't need the original, save some memory
    commonSetup = commonSetupSingleLesionLung(patient, datafile, cst, ct, pathToPredDose, "box", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;
    constraint_cst = commonSetup.constraint_cst;

    %% update CST here to remove OARs which are not close to failing constraint
    oars = {"Heart", "Esophagus", "Trachea", "BronchialTree", "SpinalCanal", "GreatVes", "Chestwall"};
    dose = datafile(strcmp(datafile.Patient, patient), :).Dose;
    dose_to_optimize_to = load(pathToPredDose).array;

    for i = 1:length(oars)
        o = oars{i};
        maxDose = getMaxDose(cst, o, dose_to_optimize_to) * pln.numOfFractions;
        constraint = getMaxConstraint(o, pln.numOfFractions, dose);        
        if maxDose < 1 * constraint
            cst{strcmp(cst(:, 2), o), 6} = {};
        end
        oar_indices_size = size(cst{strcmp(cst(:, 2), o), 4}{1});
        if oar_indices_size(1) == 0
            cst{strcmp(cst(:, 2), o), 6} = {};
        end
    end

    %Set all priorities to 1
    for i = 1:size(cst, 1)
        cst{i, 5}.Priority = 1;
    end

    cst{strcmp(cst(:, 2), 'Lung_Eval'), 6} = {}; %Don't need this for single lesion lung

    %% 

    pln.machine         = 'TBFFF_CustomFinal';

    pln2 = pln;
    pln2.propOpt.VMAToptions.maxGantryAngleSpacing    = -pln.propOpt.VMAToptions.maxGantryAngleSpacing; % Max gantry angle spacing for dose calculation
    pln2.propOpt.VMAToptions.maxDAOGantryAngleSpacing = -pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing;      % Max gantry angle spacing for DAO
    pln2.propOpt.VMAToptions.maxFMOGantryAngleSpacing = -pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing;     % Max gantry angle spacing for FMO
    pln2.propOpt.VMAToptions.startingAngle = pln.propOpt.VMAToptions.finishingAngle;
    pln2.propOpt.VMAToptions.finishingAngle = pln.propOpt.VMAToptions.startingAngle; 

    pln2 = matRad_VMATGantryAngles(pln2, cst, ct);
    stf2 = matRad_generateStf(ct,cst,pln2);

    if ~pln.propOpt.VMAToptions.continuousAperture %Don't need to do this for continuous aperture     
        stf2(1).propVMAT.timeFac = [0.5,0.5];
    end

    dij = load(fullfile(dijDir, patient)).dij;

    % Skipping calculation of dij2
    new_col_indices = [];
    num_beams = dij.numOfBeams; 
    num_rays = dij.numOfRaysPerBeam(1); 
    for beam = num_beams:-1:1 % Loop through beams in reverse order
        for ray = 1:num_rays % Loop through rays
            new_col_indices = [new_col_indices, (beam - 1) * num_rays + ray];
        end
    end
    dij2 = dij;
    dij2.physicalDose = {dij.physicalDose{1}(:, new_col_indices)};

    resultGUI1 = matRad_fluenceOptimization(dij,cst,pln,stf, 2, pathToPredDose);
    
    %% Instead of separate FMO for resultGUI2, just duplicate resultGUI1
    new_w = zeros(size(resultGUI1.w));
    new_zl = zeros(size(resultGUI1.w));
    new_zu = zeros(size(resultGUI1.w));
    
    idx = 1;
    interval = stf(1).numOfRays;
    for i = dij.totalNumOfBixels:-interval:1
        new_w(i-interval+1:i) = resultGUI1.w(idx:idx+interval-1);
        new_zl(i-interval+1:i) = resultGUI1.usedOptimizer.resultInfo.zl(idx:idx+interval-1);
        new_zu(i-interval+1:i) = resultGUI1.usedOptimizer.resultInfo.zu(idx:idx+interval-1);
        idx = idx + interval;
    end
    
    resultGUI2.w = new_w;
    resultGUI2.wUnsequenced = new_w;
    resultGUI2.usedOptimizer = resultGUI1.usedOptimizer;
    resultGUI2.usedOptimizer.wResult = new_w;
    resultGUI2.usedOptimizer.resultInfo.x = new_w;
    resultGUI2.usedOptimizer.resultInfo.zl = new_zl;
    resultGUI2.usedOptimizer.resultInfo.zu = new_zu;
    resultGUI2.info = resultGUI2.usedOptimizer.resultInfo;
    
    % Deleting physical doses to save memory
    if isfield(resultGUI1, 'physicalDose')
        resultGUI1 = rmfield(resultGUI1, 'physicalDose');
    end
    
    for i = 1:dij.numOfBeams
        fieldName = ['physicalDose_beam' num2str(i)];
        if isfield(resultGUI1, fieldName)
            resultGUI1 = rmfield(resultGUI1, fieldName);
        end
    end
    %resultGUI2 = matRad_fluenceOptimization(dij2,cst,pln2,stf2, 2, pathToPredDose);

    resultGUI1 = matRad_siochiLeafSequencing(resultGUI1,stf,dij,pln,0);
    resultGUI2 = matRad_siochiLeafSequencing(resultGUI2,stf2,dij2,pln2,0);

    resultGUI3.w = [resultGUI1.w;resultGUI2.w];
    resultGUI3.wUnsequenced = [resultGUI1.wUnsequenced;resultGUI2.wUnsequenced];
    resultGUI3.wSequenced = [resultGUI1.wSequenced;resultGUI2.wSequenced];
    
    %Create new ApertureInfo & Update vars
    newApertureInfo = resultGUI1.apertureInfo;
    newApertureInfo.jacobiScale = [resultGUI1.apertureInfo.jacobiScale;resultGUI2.apertureInfo.jacobiScale];
    newApertureInfo.totalNumOfBixels = resultGUI1.apertureInfo.totalNumOfBixels + resultGUI2.apertureInfo.totalNumOfBixels;
    newApertureInfo.totalNumOfShapes = resultGUI1.apertureInfo.totalNumOfShapes + resultGUI2.apertureInfo.totalNumOfShapes;
    newApertureInfo.totalNumOfOptBixels = resultGUI1.apertureInfo.totalNumOfOptBixels + resultGUI2.apertureInfo.totalNumOfOptBixels;
    newApertureInfo.doseTotalNumOfLeafPairs = resultGUI1.apertureInfo.doseTotalNumOfLeafPairs + resultGUI2.apertureInfo.doseTotalNumOfLeafPairs;
    newApertureInfo.totalNumOfLeafPairs = resultGUI1.apertureInfo.totalNumOfLeafPairs + resultGUI2.apertureInfo.totalNumOfLeafPairs;
    newApertureInfo.bixelWeights = [resultGUI1.apertureInfo.bixelWeights; resultGUI2.apertureInfo.bixelWeights];
    newApertureInfo.bixelJApVec = [resultGUI1.apertureInfo.bixelJApVec resultGUI2.apertureInfo.bixelJApVec];
    
    %Update ApertureInfo beams
    newBeam = resultGUI1.apertureInfo.beam;
    for i = 1:dij.numOfBeams
        newBeam(i+dij.numOfBeams) = resultGUI2.apertureInfo.beam(i);
    end
    for i=1+dij.numOfBeams:2*dij.numOfBeams
        newBeam(i).bixOffset = (i-1) * newBeam(i).numOfActiveLeafPairs + 1;
    end
    newApertureInfo.beam = newBeam;
    
    % Create new propVMAT and update apertureinfo with it
    newPropVMAT = resultGUI1.apertureInfo.propVMAT;
    newPropVMATBeam = resultGUI1.apertureInfo.propVMAT.beam;
    for i = 1:dij.numOfBeams
        newPropVMATBeam(i+dij.numOfBeams) = resultGUI2.apertureInfo.propVMAT.beam(i);
    end
    for i=1+dij.numOfBeams:2*dij.numOfBeams
        newPropVMATBeam(i).lastDAOIndex = newPropVMATBeam(i).lastDAOIndex + dij.numOfBeams;
        newPropVMATBeam(i).nextDAOIndex = newPropVMATBeam(i).nextDAOIndex + dij.numOfBeams;
        newPropVMATBeam(i).DAOIndex = newPropVMATBeam(i).DAOIndex + dij.numOfBeams;
        if pln.propOpt.VMAToptions.continuousAperture
            newPropVMATBeam(i).timeFacInd = [0 i 0];
            newPropVMATBeam(i).numOfActiveLeafPairs = 2 * newPropVMATBeam(i).numOfActiveLeafPairs;
        end
    end
    newPropVMAT.beam = newPropVMATBeam;
    newPropVMAT.jacobT = eye(2 * dij.numOfBeams);
    

    % Clear resultGUI1 and resultGUI2 to save memory
    clear resultGUI1
    clear resultGUI2
    

    if pln.propOpt.VMAToptions.continuousAperture
        newPropVMAT.numLeafSpeedConstraint = 2 * newPropVMAT.numLeafSpeedConstraint;
        newPropVMAT.numLeafSpeedConstraintDAO = 2 * newPropVMAT.numLeafSpeedConstraintDAO;
    
        %Updating timeInd -> This I'm not 100% sure about
    
        % shapeInd = 0;
        % counter = 0;
        % for i = 1:numel(newApertureInfo.beam)
        %     if newPropVMAT.beam(i).DAOBeam
        %         shapeInd = shapeInd+1;
        %         newPropVMAT.beam(i).timeInd = newApertureInfo.totalNumOfShapes+newApertureInfo.totalNumOfLeafPairs*2+shapeInd;
        %     else
        %         counter = counter + 1;
        %     end
        % end
    end
    
    newApertureInfo.propVMAT = newPropVMAT;
    
    % %Fixing vector and weight offset
    for i = 1:dij.numOfBeams*2
        newApertureInfo.beam(i).shape.vectorOffset = (dij.numOfBeams * 2 + 1) + ((i-1) * newApertureInfo.beam(1).numOfActiveLeafPairs); %Assume same numOfActiveLeafPairs throughout
        newApertureInfo.beam(i).shape.weightOffset = i;
    end
    
    
    [newApertureInfoVec, newMappingMx, newLimMx] = matRad_OptimizationProblemVMAT.matRad_daoApertureInfo2Vec(newApertureInfo);
    newApertureInfo.apertureVector = newApertureInfoVec;
    newApertureInfo.mappingMx = newMappingMx;
    newApertureInfo.limMx = newLimMx;

    resultGUI3.apertureInfo = newApertureInfo;

    %% Update plan object
    pln3 = pln;
    newPlnPropStf = pln.propStf;
    newPlnPropStf.gantryAngles = [pln.propStf.gantryAngles pln2.propStf.gantryAngles];
    newPlnPropStf.DAOGantryAngles = [pln.propStf.DAOGantryAngles pln2.propStf.DAOGantryAngles];
    newPlnPropStf.FMOGantryAngles = [pln.propStf.FMOGantryAngles pln2.propStf.FMOGantryAngles];
    newPlnPropStf.numOfBeams = pln.propStf.numOfBeams + pln2.propStf.numOfBeams;
    newPlnPropStf.couchAngles = [pln.propStf.couchAngles pln2.propStf.couchAngles];
    
    pln3.propStf = newPlnPropStf;
    
    %Create new machine, update constraints, save as temp, and set to plan
    referenceMachine = load(pln.propOpt.VMAToptions.machineConstraintFile);
    machine = referenceMachine.machine;
    
    if pln.propOpt.VMAToptions.continuousAperture
        numberOfLeafPairs = 2*resultGUI3.apertureInfo.propVMAT.numLeafSpeedConstraint*resultGUI3.apertureInfo.beam(1).numOfActiveLeafPairs;
    else
        optInd = find([resultGUI3.apertureInfo.propVMAT.beam.DAOBeam]);
        numberOfLeafPairs = 2*(numel(optInd)-1)*resultGUI3.apertureInfo.beam(1).numOfActiveLeafPairs;
    end
    
    
    %% Create new DIJ3
    dij3 = dij;
    dij3.numOfBeams = 2*dij.numOfBeams;
    dij3.numOfRaysPerBeam = [dij.numOfRaysPerBeam dij2.numOfRaysPerBeam];
    dij3.totalNumOfBixels = 2*dij.totalNumOfBixels;
    dij3.totalNumOfRays = 2*dij.totalNumOfRays;
    dij3.bixelNum = [dij.bixelNum; dij2.bixelNum];
    dij3.rayNum = [dij.rayNum; dij2.rayNum];
    dij3.beamNum = [dij.beamNum; dij.numOfBeams + dij2.beamNum];
    
    
    dij3.physicalDose = {[dij.physicalDose{1} dij2.physicalDose{1}]};
    
    clear dij
    clear dij2

    %% Calculate DAO3
    resultGUI3 = matRad_directApertureOptimization(dij3,cst,resultGUI3.apertureInfo,resultGUI3,pln3, pathToPredDose);
    toc;

    %% Repeat DAO with constraints only
    
    %modify constraint_cst to only include oars that are close to
    %failing constraints
    % oars = {"Heart", "Esophagus", "Trachea", "BronchialTree", "SpinalCanal", "GreatVes", "Chestwall"};
    % dose = datafile(strcmp(datafile.Patient, patient), :).Dose;
    % for i = 1:length(oars)
    %     o = oars{i};
    %     maxDose = getMaxDose(constraint_cst, o, resultGUI3.physicalDose) * pln.numOfFractions;
    %     constraint = getMaxConstraint(o, pln.numOfFractions, dose);        
    %     if maxDose < 0.80 * constraint
    %         constraint_cst{strcmp(constraint_cst(:, 2), o), 6} = {};
    %     end
    %     oar_indices_size = size(cst{strcmp(cst(:, 2), o), 4}{1});
    %     if oar_indices_size(1) == 0
    %         constraint_cst{strcmp(constraint_cst(:, 2), o), 6} = {};
    %     end
    % end
    % 
    % constraint_cst{strcmp(constraint_cst(:, 2), 'Lung_Eval'), 6} = {};

    %resultGUI3 = matRad_directApertureOptimization(dij3,constraint_cst,resultGUI3.apertureInfo,resultGUI3,pln3, pathToPredDose);

    %% Getting data
    % prescription = datafile(strcmp(datafile.Patient, patient), :).Dose;
    % ptv_name = datafile(strcmp(datafile.Patient, patient), :).PTVs{1};
    
    %% Scaling dose by number of fractions
    physicalDose = resultGUI3.physicalDose * pln.numOfFractions;

    % Rescaling so that PTV 95 is prescription
    % d95 = getDoseAtXPercentOfVolume(cst, ptv_name, physicalDose, 95);
    %physicalDose = physicalDose / (d95 / prescription);
    
    % disp(d95)
    % r100 = getR100(cst, ptv_name, physicalDose, prescription);
    % disp(r100)
    
    %%Calculating Dose Metrics
    %doseMetricsTable = getDoseMetrics(cst, physicalDose, ct, ptv_name, prescription, result_save_path);

    save(dose_save_path, 'physicalDose');


end


function r100 = getR100(cst, ptvName, doseCube, prescription)
indices     = cst{strcmp(cst(:, 2), ptvName), 4}{1};
numVoxelsGreaterThanPrescription = sum(sum(sum(doseCube > prescription)));
r100 = numVoxelsGreaterThanPrescription/length(indices);
end

function dXPercent = getDoseAtXPercentOfVolume(cst, oarName, doseCube, threshold)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

dXPercent = prctile(doseInVOI, (100-threshold));
end

function maxDose = getMaxDose(cst, oarName, doseCube)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

maxDose = max(doseInVOI);
end

function maxDose = getMaxDoseToXCC(cst, oarName, doseCube, singleVoxelVolumeInCC, numOfCC)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);
sortedDose = sort(doseInVOI, 'descend');
numVoxelsNeeded = int32(numOfCC / singleVoxelVolumeInCC);
maxDose = sortedDose(numVoxelsNeeded);

end

function maxConstraint = getMaxConstraint(oarName, fraction, prescription)
fraction = round(fraction);
    if fraction == 1
        scDmax = 14;
        heartDmax = 22;
        esoDmax = 15.4;
        tracheaDmax = 20.2;
        pbtDmax = 20.2;
        gvDmax = 37;
        cwDmax = 30;
    elseif fraction == 3
        scDmax = 22;
        heartDmax = 30;
        esoDmax = 27;
        tracheaDmax = 30;
        pbtDmax = 30;
        gvDmax = 45;
        cwDmax = 1.1 * prescription;
     elseif fraction == 5
        scDmax = 30;
        heartDmax = 38;
        esoDmax = 35;
        tracheaDmax = 38;
        pbtDmax = 38;
        gvDmax = 53;
        cwDmax = 1.1 * prescription;
    elseif fraction == 8
        scDmax = 32;
        heartDmax = 46;
        esoDmax = 40;
        tracheaDmax = 40;
        pbtDmax = 40;
        gvDmax = 65;
        cwDmax = 1.1 * prescription;
    end
    switch oarName
        case 'Heart'
            maxConstraint = heartDmax;
        case 'Esophagus'
            maxConstraint = esoDmax;
        case 'Trachea'
            maxConstraint = tracheaDmax;
        case 'BronchialTree'
            maxConstraint = pbtDmax;
        case 'SpinalCanal'
            maxConstraint = scDmax;
        case 'GreatVes'
            maxConstraint = gvDmax;
        case 'Chestwall'
            maxConstraint = cwDmax;
    end
end







