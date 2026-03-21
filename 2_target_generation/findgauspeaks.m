function  filteredLabelledMatrix = findgauspeaks(X,Y,Z,minArea,maxID,nROIs1)
threshold = 0.5; % Minimum Z-value to consider
minProminence = 0.1; % Minimum prominence for peaks
nROIs1 = nROIs1(2:end);

% Step 1: Smooth the data (optional)
smoothedZ = Z;
% Step 2: Detect peaks
maximaMask = imregionalmax(smoothedZ); % Binary mask of local maxima
[peakX, peakY] = find(maximaMask); % Peak coordinates
peakValues = smoothedZ(maximaMask); % Peak Z-values

% Step 3: Filter peaks based on prominence
filteredPeaks = false(size(maximaMask));
for i = 1:length(peakValues)
    % Find local region around the peak
    localRegion = smoothedZ(max(1, peakX(i)-8):min(size(Z, 1), peakX(i)+8), ...
                            max(1, peakY(i)-8):min(size(Z, 2), peakY(i)+8));
    localMin = min(localRegion(:));
    prominence = peakValues(i) - localMin;
    if prominence > minProminence
        filteredPeaks(peakX(i), peakY(i)) = true;
    end
end

% Step 4: Threshold the data
thresholdedZ = smoothedZ .* (smoothedZ > threshold);

% Step 5: Assign regions to peaks using watershed
L = watershed(-smoothedZ); % Segment using watershed
L(~(smoothedZ > threshold)) = 0; % Ignore areas below the threshold
peakLabels = L(filteredPeaks); % Map peaks to regions
L(L~=0) = 1;
% Visualization
% Original surface with peaks
% subplot(1, 2, 1);
% surf(X, Y, Z, 'EdgeColor', 'none'); colormap('jet'); hold on;
% scatter3(X(filteredPeaks), Y(filteredPeaks), smoothedZ(filteredPeaks), ...
%          50, 'r', 'filled'); % Mark peaks
% title('Original Surface with Peaks');
% xlabel('X'); ylabel('Y'); zlabel('Z'); axis tight;

% Segmented regions
%subplot(1, 2, 2);

[labeledMatrix, numBlobs] = bwlabel(L);
blobStats = regionprops(labeledMatrix, 'Centroid', 'Area'); 
filteredLabelledMatrix = zeros(size(labeledMatrix));
filteredCentroids = [];

% Loop through each blob and apply the area filter
diff = numBlobs - length(nROIs1);
if diff>0
    c = diff-1;
else
    c = 0;
end
for i = 1:numBlobs
    if blobStats(i).Area >= minArea 
        filteredCentroids = [filteredCentroids; blobStats(i).Centroid];
        filteredLabelledMatrix(labeledMatrix == i) = maxID+i;
    end
end



% imagesc(filteredLabelledMatrix); axis equal; axis tight; colormap('parula');
% title('Segmented Regions');
% xlabel('X'); ylabel('Y');

