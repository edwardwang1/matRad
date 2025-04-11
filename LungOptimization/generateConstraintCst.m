function cst = generateConstraintCst(cst, ptvs, igtvs, prescription, fraction, ct)
ccPerVoxel = ct.resolution.x / 10 * ct.resolution.y / 10 * ct.resolution.z / 10;
if ~iscell(prescription)
    prescription = {prescription};
end
max_prescription = 0;
for i = 1:length(prescription)
    if str2double(prescription{i}) > max_prescription
        max_prescription = str2double(prescription{i});
    end
end


%Create Lung_Eval - don't need to do this if we are creating this from a
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
lung_eval_LinearIdx = cst{lung_eval_index, 4};
%Remove IGTV indices
for i = 1:length(igtvs)
    igtv_name = igtvs{i}; % Get current name
    igtv_linearIdx = cst{strcmp(cst(:, 2), igtv_name), 4};
    lung_eval_LinearIdx{1}(ismember(lung_eval_LinearIdx{1}, igtv_linearIdx{1})) = [];
end
cst{lung_eval_index, 4} = lung_eval_LinearIdx;

%Intersect OARs with box
ptv_linearIdx = cst{strcmp(cst(:, 2), ptvs{1}), 4};
[ptv_x, ptv_y, ptv_z] = ind2sub(ct.cubeDim, ptv_linearIdx{1});
    
center_x = round(mean(ptv_x));
center_y = round(mean(ptv_y));
center_z = round(mean(ptv_z));

% Expand 50 mm in each direction, converted to voxels
expand_mm = 50;
expand_voxels_x = round(expand_mm / ct.resolution.x);
expand_voxels_y = round(expand_mm / ct.resolution.y);
expand_voxels_z = round(expand_mm / ct.resolution.z);

% Compute index ranges, clamping to CT bounds
x_range = max(1, center_x - expand_voxels_x) : min(ct.cubeDim(1), center_x + expand_voxels_x);
y_range = max(1, center_y - expand_voxels_y) : min(ct.cubeDim(2), center_y + expand_voxels_y);
z_range = max(1, center_z - expand_voxels_z) : min(ct.cubeDim(3), center_z + expand_voxels_z);

% Build 3D grid and get linear indices
[X, Y, Z] = ndgrid(x_range, y_range, z_range);
box_linear_idx = sub2ind(ct.cubeDim, X(:), Y(:), Z(:));

oars = {"Heart", "Esophagus", "Trachea", "BronchialTree", "SpinalCanal", "GreatVes", "Chestwall"};
for i = 1:length(oars)
    o = oars{i};
    oar_indices     = cst{strcmp(cst(:, 2), o), 4};
    new_oar_indices = {intersect(oar_indices{1}, box_linear_idx)};
    cst{strcmp(cst(:, 2), o), 4} = new_oar_indices;
end

%Adding Objectves

if fraction == 1
    lung1500cc = 7;
    lung1000cc = 7.4;
    sc1p2cc = 7;
    sc0p35cc = 10;
    scDmax = 14;
    heart15cc = 16;
    heartDmax = 22;
    eso5cc = 11.9;
    esoDmax = 15.4;
    trachea4cc = 10.5;
    tracheaDmax = 20.2;
    pbt4cc = 10.5;
    pbtDmax = 20.2;
    gv10cc = 31;
    gvDmax = 37;
    cw1cc = 22;
    cwDmax = 30;
elseif fraction == 3
    lung1500cc = 10.5;
    lung1000cc = 11.5;
    sc1p2cc = 11.1;
    sc0p35cc = 18;
    scDmax = 22;
    heart15cc = 24;
    heartDmax = 30;
    eso5cc = 21;
    esoDmax = 27;
    trachea4cc = 15;
    tracheaDmax = 30;
    pbt4cc = 15;
    pbtDmax = 30;
    gv10cc = 39;
    gvDmax = 45;
    cw1cc = 1.1 * max_prescription;
    cwDmax = 1.1 * max_prescription;
 elseif fraction == 5
    lung1500cc = 12.5;
    lung1000cc = 13.5;
    sc1p2cc = 13.5;
    sc0p35cc = 22.5;
    scDmax = 30;
    heart15cc = 32;
    heartDmax = 38;
    eso5cc = 27.5;
    esoDmax = 35;
    trachea4cc = 18;
    tracheaDmax = 38;
    pbt4cc = 18;
    pbtDmax = 38;
    gv10cc = 47;
    gvDmax = 53;
    cw1cc = 1.1 * max_prescription;
    cwDmax = 1.1 * max_prescription;
elseif fraction == 8
    lung1500cc = 12.5;
    lung1000cc = 13.5;
    sc1p2cc = 16;
    sc0p35cc = 27;
    scDmax = 32;
    heart15cc = 39;
    heartDmax = 46;
    eso5cc = 33;
    esoDmax = 40;
    trachea4cc = 21.5;
    tracheaDmax = 40;
    pbt4cc = 21.5;
    pbtDmax = 40;
    gv10cc = 58;
    gvDmax = 65;
    cwDmax = 1.1 * max_prescription;
    cw1cc = 1.1 * max_prescription;
end

lung_volume = getVolumeOfStruct(cst, "Lung_Eval", ccPerVoxel);
sc_volume = getVolumeOfStruct(cst, "SpinalCanal", ccPerVoxel);
heart_volume = getVolumeOfStruct(cst, "Heart", ccPerVoxel);
eso_volume = getVolumeOfStruct(cst, "Esophagus", ccPerVoxel);
trachea_volume = getVolumeOfStruct(cst, "Trachea", ccPerVoxel);
pbt_volume = getVolumeOfStruct(cst, "BronchialTree", ccPerVoxel);
gv_volume = getVolumeOfStruct(cst, "GreatVes", ccPerVoxel);
cw_volume = getVolumeOfStruct(cst, "Chestwall", ccPerVoxel);
ptv_vol = 0;
for i = 1:length(ptvs)
    ptv_name = ptvs{i};
    ptv_vol = ptv_vol + getVolumeOfStruct(cst, ptv_name, ccPerVoxel);
end


if ptv_vol < 20
    lungv20 = 5;
elseif ptv_vol < 40
    lungv20 = 6;
else
    lungv20 = 10;
end


% if fraction == 1
%     cwDmax = 30;
% else
%     cwDmax = 1.1 * max_prescription;
% end

s.VOIs = {'Lung_Eval', ...
    'SpinalCanal', ...
    'Heart',  ...
    'Esophagus',  ...
    'Trachea',  ...
    'BronchialTree',  ...
    'GreatVes',  ...
    'Chestwall', ...
    };

s.Parameters = {[20, lungv20]...
    [scDmax], ...
    [heartDmax], ...
    [esoDmax], ...
    [tracheaDmax], ...
    [pbtDmax], ...
    [gvDmax], ...
    [cwDmax], ...
    }; 
s.classNames = {'DoseObjectives.matRad_MaxDVH', ... %Lung
'DoseObjectives.matRad_SquaredOverdosing', ... %SC
'DoseObjectives.matRad_SquaredOverdosing', ... %Heart
'DoseObjectives.matRad_SquaredOverdosing', ... %Eso
'DoseObjectives.matRad_SquaredOverdosing', ... %Trachea
'DoseObjectives.matRad_SquaredOverdosing', ...% PBT
'DoseObjectives.matRad_SquaredOverdosing', ... %GV
'DoseObjectives.matRad_SquaredOverdosing', ... %CW
};
s.penalties = {200, ...%Lung
    200, ... %SC
    200, ... %Heart
    200, ... %Eso
    200, ... %Trachea
    200, ... %PBT
    200, ... %GV
    200, ... %CW
    };

% for i = 1:length(ptvs)
%     ptv_name = ptvs{i};
%     igtv_name = igtvs{i};
%     if isstring(prescription{i})
%         dose = str2double();
%     else
%         dose = prescription{i};
%     end
%     s.VOIs{end+1} = ptv_name;
%     s.VOIs{end+1} = igtv_name;
%     s.VOIs{end+1} = ptv_name;
%     s.Parameters{end+1} = [dose 95];
%     s.Parameters{end+1} = [dose * 1.2 0];
%     s.Parameters{end+1} = [dose * 1.5];
%     s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
%     s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
%     s.classNames{end+1} = 'DoseObjectives.matRad_SquaredOverdosing';
%     s.penalties{end+1} = 200;
%     s.penalties{end+1} = 200;
%     s.penalties{end+1} = 200;
% end

for i = 1:length(ptvs)
    ptv_name = ptvs{i};
    igtv_name = igtvs{i};
    if isstring(prescription{i})
        dose = str2double();
    else
        dose = prescription{i};
    end
    s.VOIs{end+1} = ptv_name;
    s.VOIs{end+1} = ptv_name;
    s.VOIs{end+1} = ptv_name;
    s.Parameters{end+1} = [dose 95];
    s.Parameters{end+1} = [dose * 0.9 99];
    s.Parameters{end+1} = [dose * 1.2];
    s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
    s.classNames{end+1} = 'DoseObjectives.matRad_MinDVH';
    s.classNames{end+1} = 'DoseObjectives.matRad_SquaredOverdosing';
    s.penalties{end+1} = 200;
    s.penalties{end+1} = 200;
    s.penalties{end+1} = 200;
end

for i = 1:size(cst, 1)
    cst{i, 6} = []; % Assuming cst is a cell array
end

for j = 1:size(s.VOIs, 2)
    for i = 1:size(cst, 1)
        % disp(cst{i, 2})
        % disp(s.VOIs{j})
        if strcmp(cst{i, 2}, s.VOIs{j}) % Use isequal to compare cell arrays
            tempConstraint.className = s.classNames{j}; % Use {} to access elements in cell arrays
            if size(s.Parameters{j},2) == 3
                tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2), s.Parameters{j}(3)};
            elseif size(s.Parameters{j},2) == 2
                tempConstraint.parameters = {s.Parameters{j}(1), s.Parameters{j}(2)};
            else
                tempConstraint.parameters = {s.Parameters{j}}; 
            end
            %tempConstraint.epsilon = s.epsilon{j};
            tempConstraint.penalty = s.penalties{j};
            if isempty(cst{i, 6})
                cst{i, 6} = {tempConstraint};
            elseif size(cst{i, 6}, 2) == 1
                cst{i, 6}(2) = {tempConstraint}; 
            elseif size(cst{i, 6}, 2) == 2
                cst{i, 6}(3) = {tempConstraint}; 
            end
        end
    end
end

cst{strcmp(cst(:, 2), 'GreatVes'), 3} = 'OAR';

end

function vol = getVolumeOfStruct(cst, structure, ccPerVoxel)
indices     = cst{strcmp(cst(:, 2), structure), 4}{1};
vol = length(indices) * ccPerVoxel;
end
