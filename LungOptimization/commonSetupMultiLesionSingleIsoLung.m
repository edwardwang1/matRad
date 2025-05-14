function commonSetup = commonSetupMultiLesionSingleIsoLung(patient, datafile, cst, ct, pathToDose, method, useConstraint)
    num_fractions = datafile(strcmp(datafile.Patient, patient), :).Fraction{1};
    ptv_name = datafile(strcmp(datafile.Patient, patient), :).CurrentIsoPTVs{1};
    igtv_name = datafile(strcmp(datafile.Patient, patient), :).CurrentIsoIGTVs{1};
    dose = datafile(strcmp(datafile.Patient, patient), :).Dose{1};

    num_fractions = split(num_fractions, ","); %fractions are all the same
    num_fractions = str2double(num_fractions{1});
    ptv_name = split(ptv_name, ",");
    igtv_name = split(igtv_name, ",");
    dose = split(dose, ",");

    %% Treatment Plan

    pln.radiationMode   = 'photons';   % either photons / protons / carbon
    pln.machine         = 'TBFFF_CustomFinal';
    
    pln.numOfFractions  = num_fractions;
    
    % beam geometry settings
    pln.propStf.bixelWidth = 2.5; 
    
    % optimization settings
    pln.propOpt.bioOptimization = 'none';
    pln.propOpt.runVMAT         = true;
    pln.propOpt.runDAO          = true;
    pln.propOpt.runSequencing   = true;
    pln.propOpt.preconditioner  = false; %Was true
    pln.propOpt.numLevels       = 7;
     
    pln.propOpt.VMAToptions.machineConstraintFile = [pln.radiationMode '_' pln.machine];
    
    pln.propOpt.VMAToptions.maxGantryAngleSpacing    = 4;      % Max gantry angle spacing for dose calculation
    pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 4;      % Max gantry angle spacing for DAO
    pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 28;     % Max gantry angle spacing for FMO
    
    if length(ptv_name) > 1
        [startingAngle, finishingAngle] = getArcAnglesMultiplePTVs(ct, cst, ptv_name);
    else
    [startingAngle, finishingAngle] = getArcAngles(ct, cst, ptv_name);
    end

    if finishingAngle < startingAngle
        finishingAngle = finishingAngle + 360;
    end

    % Ensure that angular range (starting angle - finishing angle) is divisible
    % by maxGantryAngleSpacing
    
    pln.propOpt.VMAToptions.startingAngle = startingAngle; 
    pln.propOpt.VMAToptions.finishingAngle = finishingAngle; 
    pln.propOpt.VMAToptions.continuousAperture = 0;
    
    pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm] 5
    pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm] 5
    pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm] 5
    
    cst = updateCST(cst, ct, ptv_name, igtv_name, dose, num_fractions, pathToDose, method, useConstraint);
    if true
        pln = setVMCParams(pln, patient);
    end
    pln = matRad_VMATGantryAngles(pln, cst, ct);
    
    %% Generate Beam Geometry STF
    stf = matRad_generateStf(ct,cst,pln);
    
    if ~pln.propOpt.VMAToptions.continuousAperture %Don't need to do this for continuous aperture     
        stf(numel(stf)).propVMAT.timeFac = [0.5,0.5];
    end

    commonSetup.cst = cst;
    commonSetup.pln = pln;
    commonSetup.stf = stf;

end



%%
function cst = updateCST(cst, ct, ptv_names, igtv_names, doses, fraction, pathToDose, method, useConstraint)
    %% Setting up CST 
    doseCube = load(pathToDose).array;
    for i = 1:size(cst, 1)
        cst{i, 6} = []; % Assuming cst is a cell array
    end

    %Get constraint CST
    %volumePerVoxel = ct.resolution.x / 10 * ct.resolution.y / 10 * ct.resolution.z / 10;
    %constraint_cst = generateConstraintCst(cst, ptv_names, igtv_names, doses, fraction, volumePerVoxel);

    %Create External_Eval
    external_eval_index = size(cst, 1) + 1;
    cst{external_eval_index, 1} = external_eval_index -1;
    cst{external_eval_index, 2} = 'External_Eval';
    cst{external_eval_index, 3} = 'OAR';
    cst{external_eval_index, 4} = cst{strcmp(cst(:, 2), 'External'), 4};
    cst{external_eval_index, 5} = cst{strcmp(cst(:, 2), 'External'), 5};
    cst{external_eval_index, 6} = [];
    external_eval_LinearIdx = cst{external_eval_index, 4};
    %Remove PTV indices
    for i = 1:length(ptv_names)
        ptv_name = ptv_names{i};
        ptv_linearIdx = cst{strcmp(cst(:, 2), ptv_name), 4};
        external_eval_LinearIdx{1}(ismember(external_eval_LinearIdx{1}, ptv_linearIdx{1})) = [];
    end

    %Update external_eval depending on method if necessary
    if strcmp(method, "naive")
    elseif strcmp(method, "box")
        %Create a 5cm box around the center of the PTV
        [x, y, z] = ind2sub(ct.cubeDim, external_eval_LinearIdx);
        center_x = round(mean(x));
        center_y = round(mean(y));
        center_z = round(mean(z));
        expand_voxels_x = round(5 / ct.resolution.x * 10); 
        expand_voxels_y = round(5 / ct.resolution.y * 10);
        expand_voxels_z = round(5 / ct.resolution.z * 10);

        x_range = max(1, center_x - expand_voxels_x) : min(array_size(1), center_x + expand_voxels_x);
        y_range = max(1, center_y - expand_voxels_y) : min(array_size(2), center_y + expand_voxels_y);
        z_range = max(1, center_z - expand_voxels_z) : min(array_size(3), center_z + expand_voxels_z);

        [X, Y, Z] = ndgrid(x_range, y_range, z_range);
        box_linear_idx = sub2ind(array_size, X(:), Y(:), Z(:));
        external_eval_LinearIdx = {intersect(external_eval_LinearIdx{1}, box_linear_idx)};

    elseif strcmp(method, "threshold")
        if fraction == 1
            threshold = 7;
        elseif fraction == 3
            threshold = 10.5;
        elseif fraction == 5
            threshold = 12.5;
        elseif fraction == 8
            threshold = 12.5;
        else
            error('Invalid value for fraction')
        end

        greater_than_threshold_linear_idx = find(doseCube >= threshold);
        external_eval_LinearIdx = {intersect(external_eval_LinearIdx{1}, greater_than_threshold_linear_idx)};


    else 
        error("Wrong Method Selected")
    end


    %This chunk of code makes external eval interesct with the isodose
    %volume of the lowest constraint
    % if fraction == 1
    %         threshold = 7;
    % elseif fraction == 3
    %     threshold = 10.5;
    % elseif fraction == 5
    %     threshold = 12.5;
    % elseif fraction == 8
    %     threshold = 12.5;
    % else
    %     error('Invalid value for fraction')
    % end
    % greater_than_threshold_linear_idx = find(doseCube >= threshold);
    % external_eval_LinearIdx = {intersect(external_eval_LinearIdx{1}, greater_than_threshold_linear_idx)};


    cst{external_eval_index, 4} = external_eval_LinearIdx;

    
    
    %Create Lung_Eval
    lung_eval_index = size(cst, 1) + 1;
    cst{lung_eval_index, 1} = lung_eval_index -1;
    cst{lung_eval_index, 2} = 'Lung_Eval';
    cst{lung_eval_index, 3} = 'OAR';
    cst{lung_eval_index, 5} = cst{strcmp(cst(:, 2), 'Lung_R'), 5};
    cst{lung_eval_index, 6} = [];
    lungR_linearIdx = cst{strcmp(cst(:, 2), 'Lung_R'), 4};
    lungL_linearIdx = cst{strcmp(cst(:, 2), 'Lung_L'), 4};
    combinedCell = [lungR_linearIdx{1}; lungL_linearIdx{1}];
    cst{lung_eval_index, 4} = {unique(combinedCell)};
    lung_eval_LinearIdx = cst{lung_eval_index, 4};
    %Remove IGTV indices
    for i = 1:length(igtv_names)
        igtv_name = igtv_names{i}; % Get current name
        igtv_linearIdx = cst{strcmp(cst(:, 2), igtv_name), 4};
        lung_eval_LinearIdx{1}(ismember(lung_eval_LinearIdx{1}, igtv_linearIdx{1})) = [];
    end
    cst{lung_eval_index, 4} = lung_eval_LinearIdx;
    
    s.VOIs = {'External_Eval'};
    s.Parameters = {20};
    s.classNames = {'DoseObjectives.matRad_SquaredDeviation'};
    s.penalties = {100};


    for i = 1:length(ptv_names)
        ptv_name = ptv_names{i};
        igtv_name = igtv_names{i};
        dose = str2double(doses{i});
        s.VOIs{end+1} = ptv_name;
        s.VOIs{end+1} = igtv_name;
        s.VOIs{end+1} = ptv_name;
        s.Parameters{end+1} = [dose 95];
        s.Parameters{end+1} = [dose * 1.2 0];
        s.Parameters{end+1} = [dose * 1.5];
        s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
        s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
        s.classNames{end+1} = 'DoseObjectives.matRad_SquaredOverdosing';
        s.penalties{end+1} = 200;
        s.penalties{end+1} = 200;
        s.penalties{end+1} = 200;
    end


    for j = 1:size(s.VOIs, 2)
        for i = 1:size(cst, 1)
            % disp(cst{i, 2})
            % disp(s.VOIs{j})
            if strcmp(cst{i, 2}, s.VOIs{j}) % Use isequal to compare cell arrays
                tempConstraint.className = s.classNames{j}; % Use {} to access elements in cell arrays
                if size(s.Parameters{j},2) == 2
                    tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2)};
                else
                    tempConstraint.parameters = {s.Parameters{j}}; 
                end
                tempConstraint.penalty = s.penalties{j};
                if isempty(cst{i, 6})
                    cst{i, 6} = {tempConstraint};
                else
                    cst{i, 6}(2) = {tempConstraint}; 
                end
            end
        end
    end

    if useConstraint
        disp("Todo")
    end

    % for i = 1:length(cst)
    %     if strcmp(cst{i,2}, 'GreatVessels')
    %         cst{i, 3} = 'OAR';
    %     elseif strcmp(cst{i,2}, 'GreatVes')
    %         cst{i, 3} = 'OAR';
    %     end
    % end

    %Change all types to OAR first, then only change current PTVs and ISOs
    %to TARGET

    for i = 1:length(cst)
        cst{i, 3} = 'OAR';
    end

    for i = 1:length(cst)
        curr_name = cst{i, 2};
        if ismember(curr_name, ptv_names)
            cst{i, 3} = 'TARGET';
        elseif ismember(curr_name, igtv_names)
            cst{i, 3} = 'TARGET';
        end
    end

end

function maxDose = getMaxDose(cst, oarName, doseCube)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

maxDose = max(doseInVOI);
end


function isFail = verifyMaxDosePass(oarName, maxDose, fraction, prescription)
isFail = -1;
if strcmp(oarName, 'SpinalCanal')
    if fraction == 1
        isFail = maxDose >= 14;
    elseif fraction == 3
        isFail = maxDose >= 22;
    elseif fraction == 5
        isFail = maxDose >= 30;
    elseif fraction == 8
        isFail = maxDose >= 32;
    end
elseif strcmp(oarName, 'Heart')
     if fraction == 1
        isFail = maxDose >= 22;
    elseif fraction == 3
        isFail = maxDose >= 30;
    elseif fraction == 5
        isFail = maxDose >= 38;
    elseif fraction == 8
        isFail = maxDose >= 46;
     end   
elseif strcmp(oarName, 'Esophagus')
     if fraction == 1
        isFail = maxDose >= 15.4;
    elseif fraction == 3
        isFail = maxDose >= 27;
    elseif fraction == 5
        isFail = maxDose >= 35;
    elseif fraction == 8
        isFail = maxDose >= 40;
    end   
elseif strcmp(oarName, 'Trachea')
     if fraction == 1
        isFail = maxDose >= 20.2;
    elseif fraction == 3
        isFail = maxDose >= 30;
    elseif fraction == 5
        isFail = maxDose >= 38;
    elseif fraction == 8
        isFail = maxDose >= 40;
    end   
elseif strcmp(oarName, 'GreatVes')
     if fraction == 1
        isFail = maxDose >= 37;
    elseif fraction == 3
        isFail = maxDose >= 45;
    elseif fraction == 5
        isFail = maxDose >= 53;
    elseif fraction == 8
        isFail = maxDose >= 65;
     end   
elseif strcmp(oarName, 'BronchialTree')
     if fraction == 1
        isFail = maxDose >= 20.2;
    elseif fraction == 3
        isFail = maxDose >= 30;
    elseif fraction == 5
        isFail = maxDose >= 38;
    elseif fraction == 8
        isFail = maxDose >= 40;
    end   
 elseif strcmp(oarName, 'Chestwall')
     if fraction == 1
        isFail = maxDose >= 30;
    elseif fraction == 3
        isFail = false;
    elseif fraction == 5
        isFail = false;
    elseif fraction == 8
        isFail = false;
     end   
end
 if isFail == -1
     error('Invalid combo of oar and fraction')
end
 

end


function pln = setVMCParams(pln, patient)
    pln.propDoseCal = 'vmc';
    pln.propDoseCalc.vmc = 1;
    pln.propDoseCalc.vmcOptions.source = 'phsp';
    pln.propDoseCalc.vmcOptions.version = 'Carleton';
    pln.propDoseCalc.vmcOptions.phspBaseName = append(patient, '_Varian_6FFF_F5-EW');
    pln.propDoseCalc.vmcOptions.SCD = 550; %550 is our truebeam data
    pln.propDoseCalc.vmcOptions.SAD = 1000;
    pln.propDoseCalc.vmcOptions.dumpDose = 1;
end
