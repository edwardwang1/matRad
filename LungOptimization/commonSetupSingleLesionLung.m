function commonSetup = commonSetupSingleLesionLung(patient, datafile, cst, ct, pathToDose, method, useConstraint)
    num_fractions = datafile(strcmp(datafile.Patient, patient), :).Fraction;
    ptv_name = datafile(strcmp(datafile.Patient, patient), :).PTVs{1};
    igtv_name = datafile(strcmp(datafile.Patient, patient), :).IGTVs{1};
    dose = datafile(strcmp(datafile.Patient, patient), :).Dose;

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
    
    
    [startingAngle, finishingAngle] = getArcAngles(ct, cst, ptv_name);

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
    pln.propDoseCalc.doseGrid.resolution.z = 3  ; % [mm] 5
    
    [cst, constraint_cst] = updateCST(cst, ct, ptv_name, igtv_name, dose, num_fractions, pathToDose, method, useConstraint);
    pln = matRad_VMATGantryAngles(pln, cst, ct);
    
    %% Generate Beam Geometry STF
    stf = matRad_generateStf(ct,cst,pln);
    
    if ~pln.propOpt.VMAToptions.continuousAperture %Don't need to do this for continuous aperture     
        stf(numel(stf)).propVMAT.timeFac = [0.5,0.5];
    end

    commonSetup.cst = cst;
    commonSetup.pln = pln;
    commonSetup.stf = stf;
    commonSetup.constraint_cst = constraint_cst;

end



%%
function [cst, constraint_cst]  = updateCST(cst, ct, ptv_name, igtv_name, dose, fraction, pathToDose, method, useConstraint)
    %% Setting up CST 
    doseCube = load(pathToDose).array;
    for i = 1:size(cst, 1)
        cst{i, 6} = []; % Assuming cst is a cell array
    end

    if ~iscell(ptv_name)
        ptv_name = {ptv_name};
    end

    if ~iscell(igtv_name)
        igtv_name = {igtv_name};
    end

    %Get constraint CST
    [constraint_cst, constraint_s] = generateConstraintCst(cst, ptv_name, igtv_name, dose, fraction, ct);

    %Create External_Eval
    external_eval_index = size(cst, 1) + 1;
    cst{external_eval_index, 1} = external_eval_index -1;
    cst{external_eval_index, 2} = 'External_Eval';
    cst{external_eval_index, 3} = 'OAR';
    cst{external_eval_index, 4} = cst{strcmp(cst(:, 2), 'External'), 4};
    cst{external_eval_index, 5} = cst{strcmp(cst(:, 2), 'External'), 5};
    cst{external_eval_index, 6} = [];

    ptv_linearIdx = cst{strcmp(cst(:, 2), ptv_name), 4};
    if ~strcmp(method, "naive")
        %Remove PTV indices
        external_eval_LinearIdx = cst{external_eval_index, 4};
        external_eval_LinearIdx{1}(ismember(external_eval_LinearIdx{1}, ptv_linearIdx{1})) = []; 
    end

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
    %Remove IGTV indices
    igtv_linearIdx = cst{strcmp(cst(:, 2), igtv_name), 4};
    lung_eval_LinearIdx = cst{lung_eval_index, 4};
    lung_eval_LinearIdx{1}(ismember(lung_eval_LinearIdx{1}, igtv_linearIdx{1})) = [];
    cst{lung_eval_index, 4} = lung_eval_LinearIdx;

    %Create 1cm ring
    ring1cm_index = size(cst, 1) + 1;
    cst{ring1cm_index, 1} = ring1cm_index -1;
    cst{ring1cm_index, 2} = 'Ring1cm';
    cst{ring1cm_index, 3} = 'OAR';
    cst{ring1cm_index, 5} = cst{strcmp(cst(:, 2), 'External'), 5};
    % cst{ring1cm_index, 5}.Priority = 3;     %Change priorty to 3
    cst{ring1cm_index, 6} = [];
    ring1cm_LinearIdx = {getExpansionMM(ct, 10, ptv_linearIdx{1})};
    ring1cm_LinearIdx{1}(ismember(ring1cm_LinearIdx{1}, ptv_linearIdx{1})) = []; 
    cst{ring1cm_index, 4} = ring1cm_LinearIdx;
    
    %Update external_eval depending on method if necessary
    if strcmp(method, "naive") || strcmp(method, "hybrid")
        external_eval_LinearIdx = cst{external_eval_index, 4};
    elseif strcmp(method, "box")
        %Create a 50mm box around the center of the PTV
        [ptv_x, ptv_y, ptv_z] = ind2sub(ct.cubeDim, ptv_linearIdx{1});
    
        center_x = round(mean(ptv_x));
        center_y = round(mean(ptv_y));
        center_z = round(mean(ptv_z));

        % Expand 50 mm in each direction, converted to voxels
        expand_mm = 50;
        expand_voxels_x = round(expand_mm / ct.resolution.x);
        expand_voxels_y = round(expand_mm / ct.resolution.y);
        expand_voxels_z = round(expand_mm / ct.resolution.z);
        
        % Compute index ranges, clamping to CT bounds
        x_range = max(1, center_x - expand_voxels_x) : min(ct.cubeDim(1), center_x + expand_voxels_x);
        y_range = max(1, center_y - expand_voxels_y) : min(ct.cubeDim(2), center_y + expand_voxels_y);
        z_range = max(1, center_z - expand_voxels_z) : min(ct.cubeDim(3), center_z + expand_voxels_z);
        
        % Build 3D grid and get linear indices
        [X, Y, Z] = ndgrid(x_range, y_range, z_range);
        box_linear_idx = sub2ind(ct.cubeDim, X(:), Y(:), Z(:));
        external_eval_LinearIdx = {intersect(external_eval_LinearIdx{1}, box_linear_idx)};

        %Also intersect all the OARs with the box
        oars = {"Heart", "Esophagus", "Trachea", "BronchialTree", "SpinalCanal", "GreatVes", "Chestwall"};
        for i = 1:length(oars)
            o = oars{i};
            oar_indices     = cst{strcmp(cst(:, 2), o), 4};
            new_oar_indices = {intersect(oar_indices{1}, box_linear_idx)};
            cst{strcmp(cst(:, 2), o), 4} = new_oar_indices;
        end
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

    elseif strcmp(method, "overwrite")
        %Overwrite OARs that fail (except lungs) with their hard
        %constraint
        %check every oar except_lung
        oars_to_check = {'SpinalCanal', 'Heart', 'Esophagus', 'Trachea', 'BronchialTree', 'GreatVes', 'Chestwall'};

        for i = 1:length(oars_to_check)
            current_oar = oars_to_check{i};
            dMax = getMaxDose(cst, current_oar, doseCube);
            %disp(current_oar)
            %disp(dMax)
            if verifyMaxDosePass(current_oar, dMax, fraction, dose) %checks if fails
                cst{strcmp(cst(:, 2), current_oar), 6} = constraint_cst{strcmp(constraint_cst(:, 2), current_oar), 6};
                current_oar_idx = cst{strcmp(cst(:, 2), current_oar), 4};
                external_eval_LinearIdx{1}(ismember(external_eval_LinearIdx{1}, current_oar_idx{1})) = [];
            end
        end


    else 
        error("Wrong Method Selected")
    end


    %This chunk of code makes external eval interesct with the isodose
    %volume of the lowest constraint


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


    cst{external_eval_index, 4} = external_eval_LinearIdx;
    
    
    % 95% of ptv receives prescription
    % s.VOIs = {'External_Eval', ptv_name, igtv_name, ptv_name};
    % s.Parameters = {20, [dose 95], [dose 100], [dose * 1.5  1]};
    % s.classNames = {'DoseObjectives.matRad_SquaredDeviation', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH'};
    % s.penalties = {100, 300, 300, 100};

    % s.VOIs = {'External_Eval', ptv_name, ptv_name, 'Ring1cm'};
    % s.Parameters = {20, [dose*1.03 95 100], [dose * 1.5  0], dose};
    % s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseConstraints.matRad_MinMaxDVH', 'DoseObjectives.matRad_MaxDVH', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {10, 100, 100, 100};

    % s.VOIs = {'External_Eval', ptv_name, ptv_name, 'Ring1cm'};
    % s.Parameters = {20, [dose*1.03 95], [dose * 1.5  0], dose};
    % s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {10, 200, 100, 100};

    % s.VOIs = {'External_Eval', ptv_name, ptv_name, 'Ring1cm'};
    % s.Parameters = {20, [dose*1.03 95], [dose * 1.5  0], dose};
    % s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {100, 3000, 2000, 500};

    if strcmp(method, "hybrid")
        s.VOIs = {'External_Eval', ptv_name, ptv_name, 'Ring1cm', 'Lung_Eval'};
        s.Parameters = {20, [dose*1.03 95], [dose * 1.5  0], dose, [20, 6]};
        s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MaxDVH'};
        s.penalties = {100, 3000, 2000, 500, 2000};
    elseif strcmp(method, "naive")
        s.VOIs = {'External_Eval', ptv_name, ptv_name};
        s.Parameters = {20, [dose*1.03 95], [dose * 1.5  0]};
        s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH'};
        s.penalties = {100, 200, 200,};
    else
        error("Not yet implemented")
    end


    %Use squared overdosing for external_eval because we don't care if
    %optimized dose is better than predicted dose
    % s.VOIs = {'External_Eval', ptv_name, igtv_name, ptv_name};
    % s.Parameters = {20, [dose 95], [dose * 1.2 0], [dose * 1.5]};
    % s.classNames = {'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {100, 200, 200, 100};
    
    % s.VOIs = {'External_Eval', ptv_name, ptv_name, 'Esophagus', 'Esophagus'};
    % s.Parameters = {20, [dose 95], dose * 1.2, 40, [0, 40, 1]};
    % s.classNames = {'DoseObjectives.matRad_SquaredDeviation', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseConstraints.matRad_MinMaxDose'};
    % s.penalties = {10, 1000, 1000, 1000, 0.001};

    % s.VOIs = {'External_Eval', ptv_name, 'Chestwall', 'Esophagus'};
    % s.Parameters = {20, [dose*1.05 95 100], dose, 0.9*40};
    % s.classNames = {'DoseObjectives.matRad_SquaredDeviation', 'DoseConstraints.matRad_MinMaxDVH', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {10, 1000, 1000, 1000};

    % s.VOIs = {ptv_name, ptv_name, igtv_name, ptv_name, 'Ring1cm'};
    % s.Parameters = {[dose 95], [dose*0.9 99], [dose * 1.2], [dose * 1.2] [dose * 1.05]};
    % s.classNames = {'DoseObjectives.matRad_MinDVH','DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_SquaredOverdosing', 'DoseObjectives.matRad_SquaredOverdosing'};
    % s.penalties = {200, 200, 200, 200, 50};


    % s.VOIs = {'External_Eval', ptv_name};
    % s.Parameters = {20, 1};
    % s.classNames = {'DoseObjectives.matRad_SquaredDeviation', 'DoseObjectives.matRad_SquaredUnderdosing'};
    % s.penalties = {100, 1};

    if strcmp(method, "hybrid")
        %Add OAR constraints based on constraint CST
        for j = 1:size(constraint_s.VOIs, 2)
            s.VOIs{end + 1} = constraint_s.VOIs{j};
            s.Parameters{end + 1} = constraint_s.Parameters{j};
            s.classNames{end + 1} = constraint_s.classNames{j};
            s.penalties{end + 1} = constraint_s.penalties{j};
        end
    end
    
    
    for j = 1:size(s.VOIs, 2)
        for i = 1:size(cst, 1)
            % disp(cst{i, 2})
            % disp(s.VOIs{j})
            if strcmp(cst{i, 2}, s.VOIs{j}) % Use isequal to compare cell arrays
                tempConstraint.className = s.classNames{j}; % Use {} to access elements in cell arrays
                if size(s.Parameters{j},2) == 3
                    tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2), s.Parameters{j}(3)};
                elseif size(s.Parameters{j},2) == 2
                    tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2)};
                else
                    tempConstraint.parameters = {s.Parameters{j}}; 
                end
                if contains(tempConstraint.className, 'Constraint') && contains(tempConstraint.className, 'DVH')
                    tempConstraint.epsilon = s.penalties{j};
                elseif contains(tempConstraint.className, 'Constraint') && contains(tempConstraint.className, 'Dose')
                    tempConstraint.voxelScalingRatio = 1;
                    tempConstraint.referenceScalingVal = 0.01;
                else
                    tempConstraint.penalty = s.penalties{j};
                end
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

    

    cst{strcmp(cst(:, 2), 'GreatVes'), 3} = 'OAR';
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

function expandedLinearIdx = getExpansionMM(ct, radius, ptv_indices)
    spacing = [ct.resolution.x, ct.resolution.y, ct.resolution.z]; % [z_spacing, y_spacing, x_spacing]

    % Calculate voxel-space radii
    radii = radius ./ spacing;
    
    % Create a 3D grid for the structuring element
    
    [x, y, z] = ndgrid(-radii(1):radii(1), -radii(2):radii(2), -radii(3):radii(3));
    
    % Ellipsoid equation to create the structuring element
    structuringElement = (x / radii(1)).^2 + (y / radii(2)).^2 + (z / radii(3)).^2 <= 1;

    indices     = ptv_indices;
    binary_array = zeros(ct.cubeDim);
    binary_array(indices) = 1;
    
    dilation = fastDilate(binary_array, structuringElement);
    expandedLinearIdx = find(dilation > 0);

end