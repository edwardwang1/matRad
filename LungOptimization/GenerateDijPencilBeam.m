%% Read csv file
tic;
datafile = readtable('SingleLesionAllEQD2.csv');
patients = readlines('test.txt');
%%
for i = 1:numel(patients)
    patient = patients(i);
    dij = setup(patient, datafile);
    savefilename = fullfile('dijPBK', patient + '.mat');
    save(savefilename, 'dij', '-v7.3');
    % break
end


%%
function dij = setup(patient, datafile)
    patientParentDir = "SingleLesionLungPatients";
    matRad_rc
    load(fullfile(patientParentDir, patient));
    commonSetup = commonSetupSingleLesionLung(patient, datafile, cst, ct);
    pln = commonSetup.pln;
    stf = commonSetup.stf;
    cst = commonSetup.cst;

    pln.machine         = 'TBFFF_CustomFinal';

    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
end

%%
