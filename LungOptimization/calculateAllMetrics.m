%% Read csv file
datafile = readtable('SingleLesionAllEQD2.csv');
patients = readlines('test.txt');
doseParentDir = "E:/matRadData/SingleLesionLungDoses";
resultsDir = "LungOptimization/Results/DoseMetricsCSVs/SingleLesionPhysicalDosesOptimizedFromRawPredictionNaive";
%doseSaveDir = "SingleLesionPhysicalDoses/OverwriteFailedOars";
doseSaveDir = "E:/matRadData/SingleLesionPhysicalDosesOptimizedFromRawPredictionNaive";
patientParentDir = "E:/matRadData/SingleLesionLungPatients";


for i = 1:numel(patients)
    close all
    patient = patients(i);
    
    matRad_rc
    load(fullfile(patientParentDir, patient));

    %CST doesn't depend on what the what the predicted dose is, only need path as dummy
    %argument
    commonSetup = commonSetupSingleLesionLung(patient, datafile, cst, ct, fullfile(doseParentDir, "GAN" + patient), "naive", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;

    if ~isfield(resultGUI, "physicalDose")
        resultGUI.physicalDose = resultGUI.RBExDose;
    end
    prescription = datafile(strcmp(datafile.Patient, patient), :).Dose;
    ptv_name = datafile(strcmp(datafile.Patient, patient), :).PTVs{1};

    groundTruthDose = resultGUI.physicalDose;
    
    %have to convert back to physicalDose
    alpha_beta_ratio = 3;
    fraction = pln.numOfFractions;
    groundTruthDose = (-alpha_beta_ratio + sqrt(alpha_beta_ratio*alpha_beta_ratio + 4/fraction *(2 * groundTruthDose + alpha_beta_ratio * groundTruthDose)))/(2/fraction);

    ganPredDose = load(fullfile(doseParentDir, "GAN" + patient)).array;
    unetPredDose = load(fullfile(doseParentDir, "Unet" + patient)).array;
    hdUnetPredDose = load(fullfile(doseParentDir, "HDUnet" + patient)).array;

    ganOptDose = load(fullfile(doseSaveDir, "GAN" + patient)).physicalDose;
    unetOptDose = load(fullfile(doseSaveDir, "Unet" + patient)).physicalDose;
    %ganOptDose = load(fullfile(doseSaveDir, "HDUnet" + patient)).physicalDose;  %For now, only look at HDUnet
    %unetOptDose = load(fullfile(doseSaveDir, "HDUnet" + patient)).physicalDose;
    hdUnetOptDose = load(fullfile(doseSaveDir, "HDUnet" + patient)).physicalDose;

    optimizedConstraintOnlyDose = load(fullfile("E:\matRadData\SingleLesionPhysicalDosesOptimizedFromConstraintsOnly", "HDUnet" + patient)).physicalDose;

    %Scale so PTV95 is equal to prescription
    d95 = getDoseAtXPercentOfVolume(cst, ptv_name, ganOptDose, 95);
    ganOptDose = ganOptDose / (d95 / prescription);

    d95 = getDoseAtXPercentOfVolume(cst, ptv_name, unetOptDose, 95);
    unetOptDose = unetOptDose / (d95 / prescription);

    d95 = getDoseAtXPercentOfVolume(cst, ptv_name, hdUnetOptDose, 95);
    hdUnetOptDose = hdUnetOptDose / (d95 / prescription);

    %Precalculate linear indices for D2cm to save time
    spacing = [ct.resolution.x/10, ct.resolution.y/10, ct.resolution.z/10]; % [z_spacing, y_spacing, x_spacing]
    real_radius = 2;          % Units are in cm
    % Calculate voxel-space radii
    radii = real_radius ./ spacing;
    % Create a 3D grid for the structuring element
    [x, y, z] = ndgrid(-radii(1):radii(1), -radii(2):radii(2), -radii(3):radii(3));
    % Ellipsoid equation to create the structuring element
    structuringElementFor2Dcm = (x / radii(1)).^2 + (y / radii(2)).^2 + (z / radii(3)).^2 <= 1;
    indices     = cst{strcmp(cst(:, 2), ptv_name), 4}{1};
    binary_array = zeros(size(groundTruthDose));
    binary_array(indices) = 1;
    
    %dilation = imdilate(binary_array, structuringElementFor2Dcm);
    dilation = fastDilate(binary_array, structuringElementFor2Dcm);
    
    linearIndicesFor2Dcm = find(dilation > 0);
    
    doseMetricsTableGroundTruth = getDoseMetrics(cst, groundTruthDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableGroundTruth.("Experiment") = "GroundTruth";
    doseMetricsTableGANPred = getDoseMetrics(cst, ganPredDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableGANPred.("Experiment") = "GanPred";
    doseMetricsTableUnetPred = getDoseMetrics(cst, unetPredDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableUnetPred.("Experiment") = "UnetPred";
    doseMetricsTableHDUnetPred = getDoseMetrics(cst, hdUnetPredDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableHDUnetPred.("Experiment") = "HDUnetPred";
    doseMetricsTableGANOpt = getDoseMetrics(cst, ganOptDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableGANOpt.("Experiment") = "GanOpt";
    doseMetricsTableUnetOpt = getDoseMetrics(cst, unetOptDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableUnetOpt.("Experiment") = "UnetOpt";
    doseMetricsTableHDUnetOpt = getDoseMetrics(cst, hdUnetOptDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableHDUnetOpt.("Experiment") = "HDUnetOpt";
    doseMetricsTableConstraintOnlyOpt = getDoseMetrics(cst, optimizedConstraintOnlyDose, ct, ptv_name, prescription, "temp.csv", linearIndicesFor2Dcm);
    doseMetricsTableConstraintOnlyOpt.("Experiment") = "OptFromConstraint";
    
    results_save_path = fullfile(resultsDir, patient + ".csv"); 
    combinedTable = [doseMetricsTableGroundTruth; doseMetricsTableGANPred; doseMetricsTableUnetPred;  doseMetricsTableHDUnetPred; doseMetricsTableGANOpt; doseMetricsTableUnetOpt; doseMetricsTableHDUnetOpt; doseMetricsTableConstraintOnlyOpt];

    writetable(combinedTable, results_save_path);

end

function dXPercent = getDoseAtXPercentOfVolume(cst, oarName, doseCube, threshold)
indices     = cst{strcmp(cst(:, 2), oarName), 4}{1};
doseInVOI = doseCube(indices);

dXPercent = prctile(doseInVOI, (100-threshold));
end