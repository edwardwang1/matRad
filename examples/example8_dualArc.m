%% Example Photon Treatment Plan with VMAT direct aperture optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to input necessary parameters in the pln structure
% (iii) how to setup a photon dose calculation
% (iv) how to inversely optimize fluence directly from command window in MatLab.
% (v) how to apply a sequencing algorithm
% (vi) how to run a VMAT direct aperture optimization
% (vii) how to visually and quantitatively evaluate the result

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the head &
% neck patient into your workspace.
matRad_rc

load('HEAD_AND_NECK.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.


matRad_rc

load('HEAD_AND_NECK.mat');

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth = 5;

% optimization settings
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runVMAT         = true;
pln.propOpt.runDAO          = true;
pln.propOpt.runSequencing   = true;
pln.propOpt.preconditioner  = true;
pln.propOpt.numLevels       = 7;
 
pln.propOpt.VMAToptions.machineConstraintFile = [pln.radiationMode '_' pln.machine];

pln.propOpt.VMAToptions.maxGantryAngleSpacing    = 4;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 4;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 28;     % Max gantry angle spacing for FMO

pln.propOpt.VMAToptions.startingAngle = 0; %
pln.propOpt.VMAToptions.finishingAngle = 84; %180
pln.propOpt.VMAToptions.continuousAperture = 0;

pln.propDoseCalc.doseGrid.resolution.x = 30; % [mm] 5
pln.propDoseCalc.doseGrid.resolution.y = 30; % [mm] 5
pln.propDoseCalc.doseGrid.resolution.z = 30; % [mm] 5

%
pln2 = pln;
pln2.propOpt.VMAToptions.maxGantryAngleSpacing    = -4;      % Max gantry angle spacing for dose calculation
pln2.propOpt.VMAToptions.maxDAOGantryAngleSpacing = -4;      % Max gantry angle spacing for DAO
pln2.propOpt.VMAToptions.maxFMOGantryAngleSpacing = -28;     % Max gantry angle spacing for FMO
pln2.propOpt.VMAToptions.startingAngle = 84;
pln2.propOpt.VMAToptions.finishingAngle = 0; 

%%
% Generate dose calculation, DAO, and FMO angles from the parameters input
% above. FMO is performed only on the initGantryAngles set. In the DAO
% step, weights and leaf positions are optimized at the angles in the
% optGantryAngles set. Weights and leaf positions are interpolated at the
% angles in the gantryAngles set to increase the accuracy of the dose
% calculation (each iteration).

% FMO: optimize fluence on coarse subset of gantry angles
% Sequencing: select subset of apertures and spread to finer angles
% DAO: constrain for leaf speed, gantry rotation speed and MU rate

%pln = matRad_VMATGantryAnglesDualArc(pln,cst,ct);
pln = matRad_VMATGantryAngles(pln, cst, ct);
pln2 = matRad_VMATGantryAngles(pln2, cst, ct);

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
stf2 = matRad_generateStf(ct,cst,pln2);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);


%% Make a copy of dij for the reverse arc
dij2 = dij;

tempPhysicalDose = dij.physicalDose(1,1); 
physicalDose = tempPhysicalDose{1,1};
tempPhysicalDose = zeros(size(physicalDose));

for i=1:dij.numOfBeams
    correspondingBeam = dij.numOfBeams + 1 - i;
    tempPhysicalDose(:, (i-1) * dij.numOfRaysPerBeam(i) + 1 : i * dij.numOfRaysPerBeam(i)) = physicalDose(:, (correspondingBeam-1) * dij.numOfRaysPerBeam(correspondingBeam) + 1 : correspondingBeam * dij.numOfRaysPerBeam(correspondingBeam));
end

dij2.physicalDose(1,1) = {tempPhysicalDose};

%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights 
% which yield the best possible dose distribution according to the 
% predefined clinical objectives and constraints underlying the radiation 
% treatment. In VMAT, FMO is done only at the angles in the
% initGantryAngles set. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
%resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf, [], 'example8TrueOptim.mat');
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);

%% Make a copy of resultGUI for reverse arc
resultGUI2 = resultGUI;
new_w = zeros(size(resultGUI2.w));
new_wUnsequenced = zeros(size(resultGUI2.wUnsequenced));

for i = 1:dij.numOfBeams
    correspondingBeam = dij.numOfBeams + 1 - i;
    resultGUI2.(['physicalDose_beam' num2str(i)]) = resultGUI.(['physicalDose_beam' num2str(correspondingBeam)]);
    new_w((i-1) * dij.numOfRaysPerBeam(i) + 1 : i * dij.numOfRaysPerBeam(i)) = resultGUI2.w((correspondingBeam-1) * dij.numOfRaysPerBeam(correspondingBeam) + 1 : correspondingBeam * dij.numOfRaysPerBeam(correspondingBeam));
    new_wUnsequenced((i-1) * dij.numOfRaysPerBeam(i) + 1 : i * dij.numOfRaysPerBeam(i)) = resultGUI2.wUnsequenced((correspondingBeam-1) * dij.numOfRaysPerBeam(correspondingBeam) + 1 : correspondingBeam * dij.numOfRaysPerBeam(correspondingBeam));
end


resultGUI2 = matRad_siochiLeafSequencing(resultGUI2,stf2,dij2,pln2,0);

%resultGUI2 = matRad_fluenceOptimization(dij2,cst,pln2,stf2); %can be skipped
%matRadGUI;

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes. The fluence map at each angle in the initGantryAngles
% set is sequenced, with the resulting apertures spread to neighbouring
% angles from the optGantryAngles set.
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

%% MLC sequence second 

resultGUI2 = matRad_siochiLeafSequencing(resultGUI2,stf2,dij2,pln2,0);
%% 

a = resultGUI.apertureInfo.apertureVector;
b = resultGUI2.apertureInfo.apertureVector;

c = a-b;
plot(c);

%% Combine resultGUIs from both sequencing steps

resultGUI3.physicalDose = resultGUI.physicalDose + resultGUI2.physicalDose;
for i = 1:dij.numOfBeams
    resultGUI3.(['physicalDose_beam' num2str(i)]) = resultGUI.(['physicalDose_beam' num2str(i)]);
    resultGUI3.(['physicalDose_beam' num2str(i + dij.numOfBeams)]) = resultGUI2.(['physicalDose_beam' num2str(i)]); 
end
resultGUI3.w = [resultGUI.w;resultGUI.w];
resultGUI3.wUnsequenced = [resultGUI.wUnsequenced;resultGUI.wUnsequenced];
resultGUI3.wSequenced = [resultGUI.wSequenced;resultGUI.wSequenced];

%Create new ApertureInfo
newApertureInfo = resultGUI.apertureInfo;

% Create new propVMAT and update
newPropVMAT = resultGUI.apertureInfo.propVMAT;
newPropVMATBeam = resultGUI.apertureInfo.propVMAT.beam;
for i = 1:dij.numOfBeams
    newPropVMATBeam(i+dij.numOfBeams) = resultGUI2.apertureInfo.propVMAT.beam(i);
end
for i=1+dij.numOfBeams:2*dij.numOfBeams
    newPropVMATBeam(i).lastDAOIndex = newPropVMATBeam(i).lastDAOIndex + dij.numOfBeams;
    newPropVMATBeam(i).nextDAOIndex = newPropVMATBeam(i).nextDAOIndex + dij.numOfBeams;
    newPropVMATBeam(i).DAOIndex = newPropVMATBeam(i).DAOIndex + dij.numOfBeams;
end
newPropVMAT.beam = newPropVMATBeam;
newPropVMAT.jacobT = eye(2 * dij.numOfBeams);
newApertureInfo.propVMAT = newPropVMAT;

%Update other variables
newApertureInfo.jacobiScale = [resultGUI.apertureInfo.jacobiScale;resultGUI2.apertureInfo.jacobiScale];
newApertureInfo.totalNumOfBixels = resultGUI.apertureInfo.totalNumOfBixels + resultGUI2.apertureInfo.totalNumOfBixels;
newApertureInfo.totalNumOfShapes = resultGUI.apertureInfo.totalNumOfShapes + resultGUI2.apertureInfo.totalNumOfShapes;
newApertureInfo.totalNumOfOptBixels = resultGUI.apertureInfo.totalNumOfOptBixels + resultGUI2.apertureInfo.totalNumOfOptBixels;
newApertureInfo.doseTotalNumOfLeafPairs = resultGUI.apertureInfo.doseTotalNumOfLeafPairs + resultGUI2.apertureInfo.doseTotalNumOfLeafPairs;
newApertureInfo.totalNumOfLeafPairs = resultGUI.apertureInfo.totalNumOfLeafPairs + resultGUI2.apertureInfo.totalNumOfLeafPairs;
newApertureInfo.bixelWeights = [resultGUI.apertureInfo.bixelWeights; resultGUI2.apertureInfo.bixelWeights];
newApertureInfo.bixelJApVec = [resultGUI.apertureInfo.bixelJApVec resultGUI2.apertureInfo.bixelJApVec];


%newApertureInfo.apertureVector = [resultGUI.apertureInfo.apertureVector; resultGUI2.apertureInfo.apertureVector];
%newApertureInfo.mappingMx = [resultGUI.apertureInfo.mappingMx; resultGUI2.apertureInfo.mappingMx];
%newApertureInfo.limMx = [resultGUI.apertureInfo.limMx; resultGUI2.apertureInfo.limMx];



%Special modification for apertrueVector, mappingMx, limMx
newApertureInfo.apertureVector = [resultGUI.apertureInfo.apertureVector(1:dij.numOfBeams, :); ...
    resultGUI2.apertureInfo.apertureVector(1:dij.numOfBeams, :); ...
    resultGUI.apertureInfo.apertureVector(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    resultGUI2.apertureInfo.apertureVector(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    resultGUI.apertureInfo.apertureVector(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :); ...
    resultGUI2.apertureInfo.apertureVector(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :)];

newApertureInfo.limMx = [resultGUI.apertureInfo.limMx(1:dij.numOfBeams, :); ...
    resultGUI2.apertureInfo.limMx(1:dij.numOfBeams, :); ...
    resultGUI.apertureInfo.limMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    resultGUI2.apertureInfo.limMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    resultGUI.apertureInfo.limMx(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :); ...
    resultGUI2.apertureInfo.limMx(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :)];


tempMappingMx = resultGUI2.apertureInfo.mappingMx;
tempMappingMx(:, 1) = tempMappingMx(:, 1) + dij.numOfBeams;
tempMappingMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, 2) = tempMappingMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, 2) + dij.numOfBeams;

newApertureInfo.mappingMx = [resultGUI.apertureInfo.mappingMx(1:dij.numOfBeams, :); ...
    tempMappingMx(1:dij.numOfBeams, :); ...
    resultGUI.apertureInfo.mappingMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    tempMappingMx(dij.numOfBeams+1:numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    resultGUI.apertureInfo.mappingMx(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :); ...
    tempMappingMx(numel(resultGUI.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :)];


%newApertureInfo.apertureVector = [resultGUI.apertureInfo.apertureVector; resultGUI2.apertureInfo.apertureVector(dij.numOfBeams+1:end, :)];
%newApertureInfo.mappingMx = [resultGUI.apertureInfo.mappingMx; resultGUI2.apertureInfo.mappingMx(dij.numOfBeams+1:end, :)];
%newApertureInfo.limMx = [resultGUI.apertureInfo.limMx; resultGUI2.apertureInfo.limMx(dij.numOfBeams+1:end, :)];


%Create new Beam and update
newBeam = resultGUI.apertureInfo.beam;
for i = 1:dij.numOfBeams
    newBeam(i+dij.numOfBeams) = resultGUI2.apertureInfo.beam(i);
end
for i=1+dij.numOfBeams:2*dij.numOfBeams
    newBeam(i).bixOffset = (i-1) * newBeam(i).numOfActiveLeafPairs + 1;
end
newApertureInfo.beam = newBeam;

resultGUI3.apertureInfo = newApertureInfo;
resultGUI3.apertureInfo = rmfield(resultGUI3.apertureInfo, 'mappingMx');


%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an optimization approach where we 
% directly optimize aperture shapes and weights at the angles in the
% optGantryAngles set.  The gantry angle speed, leaf speed, and MU rate are
% constrained by the min and max values specified by the user.
resultGUItemp = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);

%% Calculate DAO 2
resultGUI2 = matRad_directApertureOptimization(dij2,cst,resultGUI2.apertureInfo,resultGUI2,pln2);

%% Create new DIJ3
dij3 = dij;
dij3.numOfBeams = 2*dij.numOfBeams;
dij3.numOfRaysPerBeam = [dij.numOfRaysPerBeam dij2.numOfRaysPerBeam];
dij3.totalNumOfBixels = 2*dij.totalNumOfBixels;
dij3.totalNumOfRays = 2*dij.totalNumOfRays;
dij3.bixelNum = [dij.bixelNum; dij2.bixelNum];
dij3.rayNum = [dij.rayNum; dij2.rayNum];
dij3.beamNum = [dij.beamNum; dij.numOfBeams + dij2.beamNum];
dij3.physicalDose = {[dij.physicalDose{1} dij2.physicalDose{1}]};

%% Calculate DAO3


%resultGUI3.apertureInfo = rmfield(resultGUI3.apertureInfo, 'limMx');

resultGUI3 = matRad_directApertureOptimization(dij3,cst,resultGUI3.apertureInfo,resultGUI3,pln);

%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation and display of DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);
matRad_showDVH(dvh,cst,pln);

