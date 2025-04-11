
%a = imdilate(binary_array, structuringElementFor2Dcm);
%c = fastDilate2(binary_array, structuringElementFor2Dcm);

function optimized_dilation =  fastDilate(binary_array, structuring_element)
% Find the bounding box of the region of interest
    idx = find(binary_array); % Find nonzero values in binary_array
    if isempty(idx) % Handle case with no foreground pixels
        optimized_dilation = binary_array;
        return;
    end

    [x, y, z] = ind2sub(size(binary_array), idx);

    % Determine the ROI limits
    x_min = max(1, min(x) - floor(size(structuring_element, 1) / 2));
    x_max = min(size(binary_array, 1), max(x) + floor(size(structuring_element, 1) / 2));
    y_min = max(1, min(y) - floor(size(structuring_element, 2) / 2));
    y_max = min(size(binary_array, 2), max(y) + floor(size(structuring_element, 2) / 2));
    z_min = max(1, min(z) - floor(size(structuring_element, 3) / 2));
    z_max = min(size(binary_array, 3), max(z) + floor(size(structuring_element, 3) / 2));

    % Crop the ROI
    cropped_array = binary_array(x_min:x_max, y_min:y_max, z_min:z_max);

    % Pad the cropped array
    padding = floor(size(structuring_element) / 2);
    padded_cropped = padarray(cropped_array, padding, 0);

    % Perform dilation on the padded region
    dilated_padded = imdilate(padded_cropped, structuring_element);

    % Remove padding
    dilated_cropped = dilated_padded(padding(1) + 1:end - padding(1), ...
                                     padding(2) + 1:end - padding(2), ...
                                     padding(3) + 1:end - padding(3));

    % Initialize the result
    optimized_dilation = zeros(size(binary_array), 'like', binary_array);

    % Place the dilated result back into the original array
    optimized_dilation(x_min:x_max, y_min:y_max, z_min:z_max) = dilated_cropped;

end