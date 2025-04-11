function [startingAngle, finishingAngle] = getArcAngles(ct, cst, ptv_name)
    %GETARCANGLES calculates start and stop angles (assuming clockwise)
    
    % Retrieve linear indices for the PTV and External structures
    ptv_linearIdx = cst{strcmp(cst(:, 2), ptv_name), 4};
    spinalcanal_linearIdx = cst{strcmp(cst(:, 2), "SpinalCanal"), 4};
    
    % Initialize arrays to hold left-right coordinates
    leftright_coords_ptv = zeros(1, length(ptv_linearIdx));
    leftright_coords_external = zeros(1, length(spinalcanal_linearIdx));
    
    % Extract left-right coordinates for PTV
    for i = 1:length(ptv_linearIdx{1})
        [~, xCoord, ~] = ind2sub(ct.cubeDim, ptv_linearIdx{1}(i));
        leftright_coords_ptv(i) = xCoord; % Append to array xCoord
    end
    
    % Extract left-right coordinates for External
    for i = 1:length(spinalcanal_linearIdx{1})
        [~, xCoord, ~] = ind2sub(ct.cubeDim, spinalcanal_linearIdx{1}(i));
        leftright_coords_external(i) = xCoord; % Append to array
    end
    
    % Calculate middle coordinate of the External structure
    middleCoord = mean(leftright_coords_external);
    
    % Calculate the centroid of the PTV structure
    ptvCentroid = mean(leftright_coords_ptv);
    
    % Determine starting and finishing angles based on centroid position
    if ptvCentroid < middleCoord
        startingAngle = 180;
        finishingAngle = 0;
    else
        startingAngle = 0;
        finishingAngle = 180;
    end
end
