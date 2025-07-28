% load('E:\matRadData\SingleLesionLungPatients\LP386.mat');
% a = load('E:\matRadData\SingleLesionLungDoses\GANLP386.mat');
% resultGUI.physicalDose = a.array;
% prescription = 60;
% fraction = 8;
% ptv_name = 'PTV60';


%% Read csv file
datafile = readtable('SingleLesionAllEQD2.csv');
patients = readlines('test.txt');
doseParentDir = "E:\matRadData\SingleLesionLungDoses\";
doseSaveDir = "E:\matRadData\SingleLesionLungDosesPreprocessed\";

for i = 1:numel(patients)
    close all
    patient = patients(i);
    preprocess(patient, datafile, fullfile(doseParentDir, "GAN" + patient), fullfile(doseSaveDir, "GAN" + patient))
    preprocess(patient, datafile, fullfile(doseParentDir, "Unet" + patient), fullfile(doseSaveDir, "Unet" + patient))
    preprocess(patient, datafile, fullfile(doseParentDir, "HDUnet" + patient), fullfile(doseSaveDir, "HDUnet" + patient))
end

%%
function preprocess(patient, datafile, pathToPredDose, dose_save_path)
    patientParentDir = "E:\matRadData\SingleLesionLungPatients\";
    load(fullfile(patientParentDir, patient));
    fraction = datafile(strcmp(datafile.Patient, patient), :).Fraction;
    ptv_name = datafile(strcmp(datafile.Patient, patient), :).PTVs{1};
    prescription = datafile(strcmp(datafile.Patient, patient), :).Dose;
    predDose = load(pathToPredDose).array;
    %% Step 1, find all OARs that fail DMax constraint
    d95 = getDoseAtXPercentOfVolume(cst, ptv_name, predDose, 95);
    doseCube = predDose / (d95 / prescription);
    
    
    heartMaxDose = getMaxDose(cst, "Heart", doseCube);
    esoMaxDose = getMaxDose(cst, "Esophagus", doseCube);
    tracheaMaxDose = getMaxDose(cst, "Trachea", doseCube);
    pbtMaxDose = getMaxDose(cst, "BronchialTree", doseCube);
    scMaxDose = getMaxDose(cst, "SpinalCanal", doseCube);
    gvMaxDose = getMaxDose(cst, "GreatVes", doseCube);
    cwMaxDose = getMaxDose(cst, "Chestwall", doseCube);
    
    oars = {"Heart", "Esophagus", "Trachea", "BronchialTree", "SpinalCanal", "GreatVes", "Chestwall"};
    
    oars_over_constraint = {};
    constraints = {};
    maxDoses = {};
    
    for i = 1:length(oars)
        o = oars{i};
        maxDose = getMaxDose(cst, o, doseCube);
        constraint = getMaxConstraint(o, fraction, prescription);
        
        % if maxDose > constraint
        %     disp(o)
        %     oars_over_constraint{end+1} = o;
        %     constraints{end + 1} = constraint;
        %     maxDoses{end+1} = maxDose;
        % end
        if maxDose > constraint
            oar_indices     = cst{strcmp(cst(:, 2), o), 4}{1};
            ptv_indices = cst{strcmp(cst(:, 2), ptv_name), 4}{1}; 
        
            [oar_x, oar_y, oar_z] = ind2sub(ct.cubeDim, oar_indices);
            [ptv_x, ptv_y, ptv_z] = ind2sub(ct.cubeDim, ptv_indices);
            oar_centroid = [mean(oar_x), mean(oar_y), mean(oar_z)];
            ptv_centroid = [mean(ptv_x), mean(ptv_y), mean(ptv_z)];
            oar_coords = [oar_x oar_y oar_z];
            ptv_coords = [ptv_x, ptv_y, ptv_z];
        
            %Get PTV Bounds
            ptv_max_x = max(ptv_x);
            ptv_max_y = max(ptv_y);
            ptv_max_z = max(ptv_z);
        
            ptv_min_x = min(ptv_x);
            ptv_min_y = min(ptv_y);
            ptv_min_z = min(ptv_z);
        
            expansion_in_each_dir_mm = 30;
            expansion_x = round(expansion_in_each_dir_mm / ct.resolution.x);
            expansion_y = round(expansion_in_each_dir_mm / ct.resolution.y);
            expansion_z = round(expansion_in_each_dir_mm / ct.resolution.z);
        
            cube_min_x = ptv_min_x - expansion_x;
            cube_max_x = ptv_max_x + expansion_x;
            cube_min_y = ptv_min_y - expansion_y;
            cube_max_y = ptv_max_y + expansion_y;
            cube_min_z = ptv_min_z - expansion_z;
            cube_max_z = ptv_max_z + expansion_z;
        
            % Get OAR indices within subcube to save time
            inner_x = [cube_min_x:cube_max_x]';
            inner_y = [cube_min_y:cube_max_y]';
            inner_z = [cube_min_z:cube_max_z]';
        
            cube_bounds.min_x = cube_min_x;
            cube_bounds.max_x = cube_max_x;
            cube_bounds.min_y = cube_min_y;
            cube_bounds.max_y = cube_max_y;
            cube_bounds.min_z = cube_min_z;
            cube_bounds.max_z = cube_max_z;
    
            oar_magnitudes = computeOARDistanceMagnitudesAnisotropic(ct, cube_bounds, oar_coords);
            ptv_magnitudes = computeOARDistanceMagnitudesAnisotropic(ct, cube_bounds, ptv_coords);
        
            % Step 3 Apply scaling field to dose
            % disp(o)
            % disp(constraint)
            % disp(maxDose)
            oarComponent = min(exp(0.1 * oar_magnitudes), 2) * min(constraint/maxDose, 1);
            ptvComponent = exp(0.1 * -ptv_magnitudes);
            
            alpha = 0; %No PTV
            magnitudes_only = (1-alpha) * oarComponent + alpha * ptvComponent;
            magnitudes_only = min(magnitudes_only, 1);
            
            scaling_array = ones(ct.cubeDim);
            scaling_array(cube_min_x:cube_max_x, cube_min_y:cube_max_y, cube_min_z:cube_max_z) = magnitudes_only;
            
            doseCube = doseCube .* scaling_array;
        end

        array = doseCube;
        save(dose_save_path, 'array');
    end
end

%%
% resultGUI.physicalDose = doseCube;
% resultGUI.RBExDose = a.array;
% matRadGUI;

%%
% getMaxDose(cst, "Chestwall", doseCube)
% getDoseAtXPercentOfVolume(cst, ptv_name, doseCube, 95)
% getR100(cst, ptv_name, doseCube, prescription)
%% 

function maxDose = getMaxDose(cst, oarName, doseCube)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

maxDose = max(doseInVOI);
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

function oar_magnitudes = computeOARDistanceMagnitudesAnisotropic(ct, cube_bounds, oar_surface_coords)
    % Dimensions of the subcube
    sz_x = cube_bounds.max_x - cube_bounds.min_x + 1;
    sz_y = cube_bounds.max_y - cube_bounds.min_y + 1;
    sz_z = cube_bounds.max_z - cube_bounds.min_z + 1;

    % Initialize binary volume for the subcube
    binaryVolume = false(sz_x, sz_y, sz_z);

    % Convert global OAR coords to subcube-local coords
    oar_sub_x = oar_surface_coords(:,1) - cube_bounds.min_x + 1;
    oar_sub_y = oar_surface_coords(:,2) - cube_bounds.min_y + 1;
    oar_sub_z = oar_surface_coords(:,3) - cube_bounds.min_z + 1;

    % Remove any out-of-bounds points
    valid_idx = oar_sub_x >= 1 & oar_sub_x <= sz_x & ...
                oar_sub_y >= 1 & oar_sub_y <= sz_y & ...
                oar_sub_z >= 1 & oar_sub_z <= sz_z;

    % Create binary volume
    linear_idx = sub2ind([sz_x sz_y sz_z], ...
                         oar_sub_x(valid_idx), ...
                         oar_sub_y(valid_idx), ...
                         oar_sub_z(valid_idx));
    binaryVolume(linear_idx) = true;

    % Get nearest voxel index for each point in the volume
    [~, idxMap] = bwdist(binaryVolume, 'euclidean');

    % Build coordinate grids
    [X, Y, Z] = ndgrid(1:sz_x, 1:sz_y, 1:sz_z);
    current_coords = [X(:), Y(:), Z(:)];

    [tx, ty, tz] = ind2sub([sz_x sz_y sz_z], idxMap(:));
    nearest_coords = [tx(:), ty(:), tz(:)];

    % Convert both to real-world coordinates (mm)
    current_mm = [ ...
        (current_coords(:,1) - 1) * ct.resolution.x, ...
        (current_coords(:,2) - 1) * ct.resolution.y, ...
        (current_coords(:,3) - 1) * ct.resolution.z];

    nearest_mm = [ ...
        (nearest_coords(:,1) - 1) * ct.resolution.x, ...
        (nearest_coords(:,2) - 1) * ct.resolution.y, ...
        (nearest_coords(:,3) - 1) * ct.resolution.z];

    % Compute Euclidean distance in mm
    diff_mm = current_mm - nearest_mm;
    distances_mm = sqrt(sum(diff_mm.^2, 2));

    % Reshape back to 3D volume
    oar_magnitudes = reshape(distances_mm, [sz_x, sz_y, sz_z]);
end

function dXPercent = getDoseAtXPercentOfVolume(cst, oarName, doseCube, threshold)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);
dXPercent = prctile(doseInVOI, (100-threshold));
end

function r100 = getR100(cst, ptvName, doseCube, prescription)
indices     = cst{strcmp(cst(:, 2), ptvName), 4}{1};
numVoxelsGreaterThanPrescription = sum(sum(sum(doseCube > prescription)));
r100 = numVoxelsGreaterThanPrescription/length(indices);
end