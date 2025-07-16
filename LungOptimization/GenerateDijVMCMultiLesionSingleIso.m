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
    savefilename = fullfile('~/projects/def-mattonen/edwardw1/matRadData/dijVMCMultiLesion3', patient + '.mat');
    save(savefilename, 'dij', '-v7.3');
    % break
toc;
end


%%
function dij = setup(patient, datafile)
    patientParentDir = "~/projects/def-mattonen/edwardw1/matRadData/MultiLesionSingleIsoLungPatientsResampled333";
    matRad_rc
    load(fullfile(patientParentDir, patient));
    doseParentDir = "~/projects/def-mattonen/edwardw1/matRadData/MultiLesionSingleIsoLungDoses";
    pathToDose = fullfile(doseParentDir, "GAN" + patient);
    commonSetup = commonSetupMultiLesionSingleIsoLung(patient, datafile, cst, ct, pathToDose, "naive", false);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;
    pln.machine         = 'TBFFF_CustomFinal';
    dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst, 0);
end

%%  
