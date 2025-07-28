%
phasespace = [ 
"Varian_6FFF_F10_EW_photons_noZlast-10mill",
"Varian_6FFF_F10_EW_photons_noZlast-100mill",
"Varian_6FFF_F10_EW_photons_noZlast-200mill"
];

numHistories = [1000, 10000, 100000];



%%

for i = 1:length(phasespace)
    phsp = phasespace(i);
    for j = 1:length(numHistories)
        numHistory = numHistories(j);
        performExperiment(phsp, numHistory)
    end
end



%%
function performExperiment(phsp, numHistory)

pln.propDoseCalc.vmc = 1;
pln.propDoseCalc.vmcOptions.source = 'phsp';
pln.propDoseCalc.vmcOptions.version = 'vfcc';
pln.propDoseCalc.vmcOptions.phspBaseName = phsp;
pln.propDoseCalc.vmcOptions.SCD = 533;
pln.propDoseCalc.vmcOptions.SAD = 1000;
pln.propDoseCalc.vmcOptions.dumpDose = 1;
pln.propDoseCalc.vmcOptions.nCasePerBixel = numHistory;


save_prefix = phsp + "_" + numHistory;
dij_save_name = save_prefix + "Dij";
resultGUI_save_name = save_prefix + "ResultGUIl";
disp(save_prefix)
matRad_rc;
load('vmcCUBE400mm.mat')


% meta information for treatment plan
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth      = 100; % 10x10cm field size
pln.propStf.gantryAngles    = 0; % [°]
pln.propStf.couchAngles     = 0; % [°]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%%

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below


%% define steering file by hand

%
stf = matRad_generateStf(ct,cst,pln, 0);

%% dose calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,1);

%% start gui for visualization of result
w = ones(1);

%% Calculate dose
resultGUI = matRad_calcCubes(w, dij);
save(dij_save_name, 'dij', '-v7.3');
save(resultGUI_save_name, 'resultGUI', '-v7.3')

end


