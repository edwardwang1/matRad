%% Read csv file
tic;
opts = detectImportOptions('SingleIsoEQD2.csv', 'Delimiter', ',');
opts = setvaropts(opts, {'Fraction', 'Dose'}, 'Type', 'string');


datafile = readtable('SingleIsoEQD2.csv', opts);
patients = readlines('testIso.txt');
%%
for i = 1:numel(patients)
    patient = patients(i);
    dij = setup(patient, datafile);
    savefilename = fullfile('E:\AutomatedLungSBRTPlanningData\dijPBKMultiLesion', patient + '.mat');
    save(savefilename, 'dij', '-v7.3');
    % break
end


%%
function dij = setup(patient, datafile)
    patientParentDir = "MultiLesionSingleIsoLungPatients";
    matRad_rc
    load(fullfile(patientParentDir, patient));
    doseParentDir = "MultiLesionSingleIsoLungDoses";
    pathToDose = fullfile(doseParentDir, "GAN" + patient);
    commonSetup = commonSetupMultiLesionSingleIsoLung(patient, datafile, cst, ct, pathToDose, "naive", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;

    pln.machine         = 'TBFFF_CustomFinal';

    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
end

%%  
