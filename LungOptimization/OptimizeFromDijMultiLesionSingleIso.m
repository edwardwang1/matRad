%% Read csv file
opts = detectImportOptions('SingleIsoEQD2.csv', 'Delimiter', ',');
opts = setvaropts(opts, {'Fraction', 'Dose'}, 'Type', 'string');

datafile = readtable('SingleIsoEQD2.csv', opts);

patients = readlines('testIso.txt');



doseParentDir = "MultiLesionSingleIsoLungDoses";

physicalDoseSaveDir = "MultiLesionSingleIsoPhysicalDosesScaled";
EQD2DoseSaveDir = "MultiLesionSingleIsoEQD2DosesScaled";


%%
for i = 1:numel(patients)
    close all
    patient = patients(i);
    runInversePlanning(patient, datafile, fullfile(doseParentDir, "GAN_Adv" + patient), fullfile(physicalDoseSaveDir, "GAN_Adv" + patient + ".mat"), fullfile(EQD2DoseSaveDir, "GAN_Adv" + patient + ".mat"));
end

%% 
%function resultGUI3 = runInversePlanning(patient, datafile, pathToPredDose, result_save_path, dose_save_path)
function runInversePlanning(patient, datafile, pathToPredDose, physicalDoseSaveDir, eqd2SaveDir)
    patientParentDir = "MultiLesionSingleIsoLungPatients";
    dijDir = "E:\AutomatedLungSBRTPlanningData\dijPBKMultiLesion";
    
    tic;
    matRad_rc
    load(fullfile(patientParentDir, patient));
    clear resultGUI %Don't need the original, save some memory
    commonSetup = commonSetupMultiLesionSingleIsoLung(patient, datafile, cst, ct, pathToPredDose, "naive", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;


    ptv_names = datafile(strcmp(datafile.Patient, patient), :).CurrentIsoPTVs{1};
    igtv_names = datafile(strcmp(datafile.Patient, patient), :).CurrentIsoIGTVs{1};
    doses = datafile(strcmp(datafile.Patient, patient), :).Dose{1};
    ptv_names = split(ptv_names, ",");
    igtv_names = split(igtv_names, ",");
    doses = split(doses, ",");

    volumePerVoxel = ct.resolution.x / 10 * ct.resolution.y / 10 * ct.resolution.z / 10;
    %constraint_cst = generateConstraintCst(cst, ptv_names, igtv_names, doses, pln.numOfFractions, volumePerVoxel);

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
    resultGUI1 = rmfield(resultGUI1, 'physicalDose');
    % for i = 1:dij.numOfBeams
    %    resultGUI1 = rmfield(resultGUI1, ['physicalDose_beam' num2str(i)]);
    % end

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
    % referenceMachine = load(pln.propOpt.VMAToptions.machineConstraintFile);
    % machine = referenceMachine.machine;
    % 
    % if pln.propOpt.VMAToptions.continuousAperture
    %     numberOfLeafPairs = 2*resultGUI3.apertureInfo.propVMAT.numLeafSpeedConstraint*resultGUI3.apertureInfo.beam(1).numOfActiveLeafPairs;
    % else
    %     optInd = find([resultGUI3.apertureInfo.propVMAT.beam.DAOBeam]);
    %     numberOfLeafPairs = 2*(numel(optInd)-1)*resultGUI3.apertureInfo.beam(1).numOfActiveLeafPairs;
    % end
    % 
    
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
    
    %% clear variables to save memory
    clear dij
    clear dij2
    clear commonSetup
    clear new_col_indices
    clear new_w
    clear new_zl
    clear new_zu
    clear newApertureInfo
    clear newApertureInfoVec
    clear newBeam
    clear newLimMx
    clear newMappingMx
    clear newPlnPropStf
    clear newPropVMAT
    clear newPropVMATBeam
    clear ct


    %% Calculate DAO3
    resultGUI3 = matRad_directApertureOptimization(dij3,cst,resultGUI3.apertureInfo,resultGUI3,pln3, pathToPredDose);
    %resultGUI3 = matRad_directApertureOptimization(dij3,constraint_cst,resultGUI3.apertureInfo,resultGUI3,pln3, pathToPredDose);
    
    toc;
    
       
    %% Scaling dose by number of fractions
    physicalDose = resultGUI3.physicalDose * pln.numOfFractions;
    
    divideRatio = 99;
    %Rescale so that Ptv D95 is always met for every target
    for j = 1:length(ptv_names)
    
        d95 = getDoseAtXPercentOfVolume(cst, ptv_names{j}, physicalDose, 95);

        tempRatio = d95 / str2double(doses{j});

        if tempRatio < divideRatio
            divideRatio = tempRatio;
        end
    end

    %physicalDose = physicalDose / divideRatio;
    eqd2Dose = convertToEQD2(physicalDose, pln.numOfFractions);



    save(physicalDoseSaveDir, 'physicalDose');
    save(eqd2SaveDir, 'eqd2Dose');


end

function dXPercent = getDoseAtXPercentOfVolume(cst, oarName, doseCube, threshold)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

dXPercent = prctile(doseInVOI, (100-threshold));
end

function eqd2Dose = convertToEQD2(physicalDose, fractions)
eqd2Dose = physicalDose .* (physicalDose/fractions + 3)/(2 + 3);
end









