%% Read csv file
tic;
opts = detectImportOptions('SingleIsoEQD2.csv', 'Delimiter', ',');
opts = setvaropts(opts, {'Fraction', 'Dose'}, 'Type', 'string');


datafile = readtable('SingleIsoEQD2.csv', opts);
patients = readlines('testForVMC.txt');
%%
for i = 1:numel(patients)
    patient = patients(i);
    dij = setup(patient, datafile);
    savefilename = fullfile('E:\matRadData\dijVMCMultiLesion3', patient + '.mat');
    %save(savefilename, 'dij', '-v7.3');
    % break
toc;
end


%%
function dij = setup(patient, datafile)
    patientParentDir = "E:\matRadData\MultiLesionSingleIsoLungPatientsResampled333";
    matRad_rc
    load(fullfile(patientParentDir, patient));
    doseParentDir = "E:\matRadData\MultiLesionSingleIsoLungDoses";
    pathToDose = fullfile(doseParentDir, "GAN" + patient);
    commonSetup = commonSetupMultiLesionSingleIsoLung(patient, datafile, cst, ct, pathToDose, "naive", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;
    pln.machine         = 'TBFFF_CustomFinal';
    dij = 0;
    %dij = matRad_calcPhotonDose(ct,stf,pln,cst);
end

%%  
