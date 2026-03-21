function [nearest] = findNearestNonBackground(x,y,image,ndist)


[rownon4, colnon4] = find(image ~= -1);

% Initialize the minimum distance and nearest value
min_dist = inf;
nearest_val = -1;

% Loop over each pixel with non-4 value to find the nearest
for j = 1:length(rownon4)
    % Calculate Euclidean distance
    dist = sqrt((y - rownon4(j))^2 + (x - colnon4(j))^2);
    
    % Check if this is the minimum distance found so far
    if dist < min_dist
        min_dist = dist;
        nearest_val = image(rownon4(j), colnon4(j));
    end
end

% If the nearest distance is less than 10, update the value
if min_dist < ndist
    nearest = nearest_val;
else
    nearest = -1;
end

