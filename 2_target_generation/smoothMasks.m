function mask_fine = smoothMasks(mask,x,y,X_new,Y_new)

mask_fine = zeros(size(X_new));

% Process each region separately
unique_labels = unique(mask);
unique_labels(unique_labels == 0) = []; % Remove background (if any)

for i = 1:length(unique_labels)
    label = unique_labels(i);
    
    % Create binary mask for the current region
    binary_mask = (mask == label);
    
    % Interpolate mask to finer grid
    binary_mask_fine = interp2(x,y,double(binary_mask), X_new, Y_new, 'cubic');
    
    % Apply Gaussian smoothing
    sigma = 15; % Standard deviation of Gaussian filter
    smoothed_binary_mask = imgaussfilt(binary_mask_fine, sigma);
    
    % Threshold to recover region
    final_binary_mask = smoothed_binary_mask > 0.1;
    
    % Assign the label back to the fine mask
    mask_fine(final_binary_mask) = label;
end

% Display results
%imagesc(mask_fine)

