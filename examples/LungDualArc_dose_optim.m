%% Patient Data Import
% Let's begin with a clear Matlab environment and import the head &
% neck patient into your workspace.
matRad_rc
load('LP054.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 5;

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

pln.propOpt.VMAToptions.startingAngle = 0; 
pln.propOpt.VMAToptions.finishingAngle = 359; 
pln.propOpt.VMAToptions.continuousAperture = 0;

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm] 5
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm] 5
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm] 5

%
pln2 = pln;
pln2.propOpt.VMAToptions.maxGantryAngleSpacing    = -4;      % Max gantry angle spacing for dose calculation
pln2.propOpt.VMAToptions.maxDAOGantryAngleSpacing = -4;      % Max gantry angle spacing for DAO
pln2.propOpt.VMAToptions.maxFMOGantryAngleSpacing = -28;     % Max gantry angle spacing for FMO
pln2.propOpt.VMAToptions.startingAngle = 182;
pln2.propOpt.VMAToptions.finishingAngle = 30; 

%% Setting up CST SABR_SYNC

s.VOIs = {'Esophagus', 'Heart', 'Chestwall', 'GreatVes', 'SpinalCanal', 'Trachea', 'Lung_Eval'...
    'PTV_a', 'PTV_b'};
s.Parameters = {35, 38, 57, 53, 30, 40, [13 37], [35 95], [35 95] };
s.classNames = {'DoseObjectives.matRad_SquaredDeviation' ,'DoseObjectives.matRad_SquaredDeviation' ,'DoseObjectives.matRad_SquaredDeviation'...
    'DoseObjectives.matRad_SquaredDeviation','DoseObjectives.matRad_SquaredDeviation','DoseObjectives.matRad_SquaredDeviation'...
    'DoseObjectives.matRad_MaxDVH', 'DoseObjectives.matRad_MinDVH', 'DoseObjectives.matRad_MinDVH'};
s.penalties = {100, 100, 100, 100, 100, 100, 100, 1000, 1000 };


for i = 1:size(cst, 1)
    cst{i, 6} = []; % Assuming cst is a cell array
end

for i = 1:size(cst, 1)
    for j = 1:size(s.VOIs, 2) % Iterate over the size of s.VOIs
        % disp(cst{i, 2})
        % disp(s.VOIs{j})
        if strcmp(cst{i, 2}, s.VOIs{j}) % Use isequal to compare cell arrays
            tempConstraint.className = s.classNames{j}; % Use {} to access elements in cell arrays
            if size(s.Parameters{j},2) == 2
                tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2)}; % Use {} to access elements in cell arrays
            else
                tempConstraint.parameters = {s.Parameters{j}}; % Use {} to access elements in cell arrays
            end
            tempConstraint.penalty = s.penalties{j};
            cst{i, 6} = {tempConstraint};
        end
    end
end
%% Setting Up CST

%defaultConstraint.className = 'DoseObjectives.matRad_SquaredDeviation';
%defaultConstraint.className = 'DoseObjectives.matRad_SquaredOverdosing';

%defaultConstraint.paramters = {30};
% defaultConstraint.penalty = 150;
% 
% constraints = [-1,8.5,26,12,5,26,26,26,26,-1,-1,3,-1,7,20,-1,-1];
% 
% %Change cst optim from SquaredOverdosing to SquaredDeviation
% for i = 1:size(cst,1)
% 
%     if i >= 6 && i <= 9
%         defaultConstraint.className = 'DoseObjectives.matRad_SquaredUnderdosing';
%     else
%         defaultConstraint.className = 'DoseConstraints.matRad_MinMaxMeanDose';
%     end
% 
%     % if ~isempty(cst{i,6})
%     %     cst{i,6}{1,1}.className = 'DoseObjectives.matRad_SquaredDeviation';
%     % end
%     defaultConstraint.parameters = {constraints(i)};
%     cst{i, 6} = {defaultConstraint};
%     %cst{i,3} = 'OAR';
% 
% 
%     if constraints(i) == -1
%         cst{i,6} = [];
%     end
% end
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

pln = matRad_VMATGantryAngles(pln, cst, ct);
pln2 = matRad_VMATGantryAngles(pln2, cst, ct);
%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
stf2 = matRad_generateStf(ct,cst,pln2);


% VERIFY LATER
stf(numel(stf)).propVMAT.timeFac = stf(2).propVMAT.timeFac; 
stf2(1).propVMAT.timeFac = stf2(2).propVMAT.timeFac;
 
% stf(numel(stf)).propVMAT.timeFac = [0.5,0.5];
% stf2(1).propVMAT.timeFac = [0.5,0.5];

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
dij2 = matRad_calcPhotonDose(ct,stf2,pln2,cst);


%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights 
% which yield the best possible dose distribution according to the 
% predefined clinical objectives and constraints underlying the radiation 
% treatment. In VMAT, FMO is done only at the angles in the
% initGantryAngles set. Once the optimization has finished, trigger once the GUI to
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf, 'LP054NonFFF.mat');
%resultGUI2 = matRad_fluenceOptimization(dij2,cst,pln2,stf2); %can be skipped


%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes. The fluence map at each angle in the initGantryAngles
% set is sequenced, with the resulting apertures spread to neighbouring
% angles from the optGantryAngles set.


resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);
%resultGUI2 = matRad_siochiLeafSequencing(resultGUI2,stf2,dij2,pln2,0);
%% Combine resultGUIs from both sequencing steps

test = resultGUI;
test2 = resultGUI2;

resultGUI3.physicalDose = test.physicalDose + test2.physicalDose;
for i = 1:dij.numOfBeams
    resultGUI3.(['physicalDose_beam' num2str(i)]) = test.(['physicalDose_beam' num2str(i)]);
    resultGUI3.(['physicalDose_beam' num2str(i + dij.numOfBeams)]) = test2.(['physicalDose_beam' num2str(i)]); 
end
resultGUI3.w = [test.w;test2.w];
resultGUI3.wUnsequenced = [test.wUnsequenced;test2.wUnsequenced];
resultGUI3.wSequenced = [test.wSequenced;test2.wSequenced];

%Create new ApertureInfo
newApertureInfo = test.apertureInfo;

% Create new propVMAT and update
newPropVMAT = test.apertureInfo.propVMAT;
newPropVMATBeam = test.apertureInfo.propVMAT.beam;
for i = 1:dij.numOfBeams
    newPropVMATBeam(i+dij.numOfBeams) = test2.apertureInfo.propVMAT.beam(i);
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
newApertureInfo.jacobiScale = [test.apertureInfo.jacobiScale;test2.apertureInfo.jacobiScale];
newApertureInfo.totalNumOfBixels = test.apertureInfo.totalNumOfBixels + test2.apertureInfo.totalNumOfBixels;
newApertureInfo.totalNumOfShapes = test.apertureInfo.totalNumOfShapes + test2.apertureInfo.totalNumOfShapes;
newApertureInfo.totalNumOfOptBixels = test.apertureInfo.totalNumOfOptBixels + test2.apertureInfo.totalNumOfOptBixels;
newApertureInfo.doseTotalNumOfLeafPairs = test.apertureInfo.doseTotalNumOfLeafPairs + test2.apertureInfo.doseTotalNumOfLeafPairs;
newApertureInfo.totalNumOfLeafPairs = test.apertureInfo.totalNumOfLeafPairs + test2.apertureInfo.totalNumOfLeafPairs;
newApertureInfo.bixelWeights = [test.apertureInfo.bixelWeights; test2.apertureInfo.bixelWeights];
newApertureInfo.bixelJApVec = [test.apertureInfo.bixelJApVec test2.apertureInfo.bixelJApVec];


%Special modification for apertrueVector, mappingMx, limMx
newApertureInfo.apertureVector = [test.apertureInfo.apertureVector(1:dij.numOfBeams, :); ...
    test2.apertureInfo.apertureVector(1:dij.numOfBeams, :); ...
    test.apertureInfo.apertureVector(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams - resultGUI.apertureInfo.totalNumOfLeafPairs, :); ...
    test2.apertureInfo.apertureVector(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs, :); ...
    test.apertureInfo.apertureVector(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    test2.apertureInfo.apertureVector(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    test.apertureInfo.apertureVector(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :); ...
    test2.apertureInfo.apertureVector(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :)];

newApertureInfo.apertureVector(1:dij.numOfBeams *2, :) = newApertureInfo.apertureVector(1:dij.numOfBeams *2, :) /2;

newApertureInfo.limMx = [test.apertureInfo.limMx(1:dij.numOfBeams, :); ...
    test2.apertureInfo.limMx(1:dij.numOfBeams, :); ...
    test.apertureInfo.limMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams - resultGUI.apertureInfo.totalNumOfLeafPairs, :); ...
    test2.apertureInfo.limMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs, :); ...
    test.apertureInfo.limMx(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    test2.apertureInfo.limMx(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :); ...
    test.apertureInfo.limMx(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :); ...
    test2.apertureInfo.limMx(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :)];


tempMappingMx = test2.apertureInfo.mappingMx;
tempMappingMx(:, 1) = tempMappingMx(:, 1) + dij.numOfBeams;
tempMappingMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, 2) = tempMappingMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, 2) + dij.numOfBeams;


mxAdj = repmat([46 46 0 0],length(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs),1);
 a = test.apertureInfo.mappingMx(1:dij.numOfBeams, :);
 b = tempMappingMx(1:dij.numOfBeams, :);
 c = test.apertureInfo.mappingMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams - resultGUI.apertureInfo.totalNumOfLeafPairs, :);
 d = test2.apertureInfo.mappingMx(dij.numOfBeams+1:numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs, :) +mxAdj;
 e = test.apertureInfo.mappingMx(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :);
 f = test2.apertureInfo.mappingMx(1+numel(test.apertureInfo.apertureVector) - dij.numOfBeams  - resultGUI.apertureInfo.totalNumOfLeafPairs:numel(test.apertureInfo.apertureVector) - dij.numOfBeams, :) +mxAdj;
 g = test.apertureInfo.mappingMx(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :);
 h = tempMappingMx(numel(test.apertureInfo.apertureVector) - dij.numOfBeams + 1:end, :);

 newApertureInfo.mappingMx = [a;b;c;d;e;f;g;h];



%Create new Beam and update
newBeam = test.apertureInfo.beam;
for i = 1:dij.numOfBeams
    newBeam(i+dij.numOfBeams) = test2.apertureInfo.beam(i);
end
for i=1+dij.numOfBeams:2*dij.numOfBeams
    newBeam(i).bixOffset = (i-1) * newBeam(i).numOfActiveLeafPairs + 1;
end
newApertureInfo.beam = newBeam;

resultGUI3.apertureInfo = newApertureInfo;
%resultGUI3.apertureInfo = rmfield(resultGUI3.apertureInfo, 'mappingMx');


% Update plan object
pln3 = pln;
newPlnPropStf = pln.propStf;
newPlnPropStf.gantryAngles = [pln.propStf.gantryAngles pln2.propStf.gantryAngles];
newPlnPropStf.DAOGantryAngles = [pln.propStf.DAOGantryAngles pln2.propStf.DAOGantryAngles];
newPlnPropStf.FMOGantryAngles = [pln.propStf.FMOGantryAngles pln2.propStf.FMOGantryAngles];
newPlnPropStf.numOfBeams = pln.propStf.numOfBeams + pln2.propStf.numOfBeams;
newPlnPropStf.couchAngles = [pln.propStf.couchAngles pln2.propStf.couchAngles];

pln3.propStf = newPlnPropStf;

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

% %% Create new DIJ3 sparsely
% dij3 = struct();
% 
% % Double the number of beams
% dij3.numOfBeams = 2 * dij.numOfBeams;
% 
% % Double the total number of bixels
% dij3.totalNumOfBixels = 2 * dij.totalNumOfBixels;
% 
% % Double the total number of rays
% dij3.totalNumOfRays = 2 * dij.totalNumOfRays;
% 
% % Concatenate arrays vertically without converting to dense
% dij3.bixelNum = [dij.bixelNum; dij2.bixelNum];
% dij3.rayNum = [dij.rayNum; dij2.rayNum];
% 
% % Concatenate beam numbers and update accordingly
% dij3.beamNum = [dij.beamNum; (dij.numOfBeams + dij2.beamNum)];
% 
% % Concatenate physical dose cell arrays
% dij3.physicalDose = [dij.physicalDose; dij2.physicalDose];


%% Calculate DAO3
resultGUI3 = matRad_directApertureOptimization(dij3,cst,resultGUI3.apertureInfo,resultGUI3,pln3);
%% For single arc test DAO
resultGUISingleArc = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);

%% open GUI
resultGUI1 = resultGUI;
resultGUI = resultGUISingleArc;
matRadGUI

%%
y = load('LP054NonFFF.mat');
dose_org = y.nonFFF;
dose_opti = resultGUI.physicalDose;
dose_diff = abs(dose_org = dose_opti);

%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo);

%% Indicator Calculation and display of DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);
matRad_showDVH(dvh,cst,pln);

