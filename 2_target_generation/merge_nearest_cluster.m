function updatedMatrix = merge_nearest_cluster(matrix, target_label)
    % Identify locations of the target cluster
    [rows, cols] = find(matrix == target_label);
    if isempty(rows) | isempty(cols)
        updatedMatrix = matrix;
        return 
    end
    % Find all unique neighboring labels
    neighbors = unique(get_neighbors(matrix, rows, cols));
    
    % Remove target_label and 0 (if background) from neighbors
    neighbors(neighbors == target_label | neighbors == 0) = [];
    
    if isempty(neighbors)
        warning('No adjacent clusters found to merge with.');
        updatedMatrix = matrix;
        return;
    end
    
    % Compute centroids of the target cluster and neighboring clusters
    target_centroid = mean([rows, cols], 1);
    nearest_label = find_nearest_cluster(matrix, neighbors, target_centroid);
    
    % Merge target cluster with the nearest one
    matrix(matrix == target_label) = nearest_label;
    updatedMatrix = matrix;
end

function neighbor_labels = get_neighbors(matrix, rows, cols)
    % Get unique neighboring labels
    neighbor_labels = [];
    [m, n] = size(matrix);
    offsets = [-1 0; 1 0; 0 -1; 0 1]; % Up, Down, Left, Right
    
    for i = 1:length(rows)
        for j = 1:size(offsets, 1)
            r = rows(i) + offsets(j, 1);
            c = cols(i) + offsets(j, 2);
            if r >= 1 && r <= m && c >= 1 && c <= n
                neighbor_labels = [neighbor_labels, matrix(r, c)];
            end
        end
    end
    neighbor_labels = unique(neighbor_labels);
end

function nearest_label = find_nearest_cluster(matrix, neighbor_labels, target_centroid)
    % Find the closest cluster based on centroid distance
    min_dist = inf;
    nearest_label = neighbor_labels(1);
    
    for label = neighbor_labels
        [r, c] = find(matrix == label);
        centroid = mean([r, c], 1);
        dist = norm(centroid - target_centroid);
        
        if dist < min_dist
            min_dist = dist;
            nearest_label = label;
        end
    end
end
