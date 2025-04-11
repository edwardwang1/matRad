function commonSetup = CommonSetupSingleLesionLungConstraints(patient, datafile, cst, ct)
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
    
    pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm] 5
    pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm] 5
    pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm] 5
    

    cst = updateCST(cst, ptv_name, igtv_name, dose);
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
function cst = updateCST(cst, ptv_name, igtv_name, dose)
    %% Setting up CST SABR_SYNC

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
    
    % 95% of ptv receives prescription
    s.VOIs = {'External_Eval', ptv_name, igtv_name, ptv_name};
    s.Parameters = {20, [dose 95], [dose 100], [dose * 1.5  1]};
    s.classNames = {'DoseObjectives.matRad_SquaredDeviation', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MaxDVH'};
    s.penalties = {100, 300, 300, 100};


    for i = 1:size(cst, 1)
        cst{i, 6} = []; % Assuming cst is a cell array
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

    cst{strcmp(cst(:, 2), 'GreatVes'), 3} = 'OAR';
end
