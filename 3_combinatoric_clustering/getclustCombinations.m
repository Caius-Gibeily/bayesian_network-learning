function allCuts = getclustCombinations(Z, minLeaves)
  
    if nargin < 2
        minLeaves = 5;
    end

    n = size(Z, 1) + 1;      % number of leaves
    nNodes = n + size(Z, 1); % total nodes (leaves + merges)
    children = cell(1, nNodes);
    heights  = zeros(1, nNodes);

    % 
    for i = 1:size(Z, 1)
        parent = n + i;
        children{parent} = [Z(i,1), Z(i,2)];
        heights(parent) = Z(i,3);
    end
    childrenPre = children;

    % 
    leafCount = zeros(1, nNodes);
    for node = 1:nNodes
        leafCount(node) = numel(getLeaves(node));
    end

    % 
    validNodes = leafCount >= minLeaves;

    % 
    allCutsRaw = getCutsRestricted(nNodes);

    % 
    keyStrings = cellfun(@(x) sprintf('%d_', sort(cell2mat(x))), allCutsRaw, 'UniformOutput', false);
    [~, ia] = unique(keyStrings, 'stable');
    allCutsRaw = allCutsRaw(ia);

    %
    allLeaves = cell(size(allCutsRaw));
    allLeafCounts = cell(size(allCutsRaw));
    for i = 1:numel(allCutsRaw)
        combo = cell2mat(allCutsRaw{i});
        leavesPerCluster = cell(size(combo));
        for j = 1:numel(combo)
            leaves = getLeaves(combo(j));
            leavesPerCluster{j} = sort(unique(leaves));
        end
        allLeaves{i} = leavesPerCluster;
        allLeafCounts{i} = cellfun(@numel, leavesPerCluster);
    end

    % 
    allHeights = cell(size(allCutsRaw));
    for i = 1:numel(allCutsRaw)
        comboHeights = cell(1, numel(allCutsRaw{i}));
        for j = 1:numel(allCutsRaw{i})
            ZID = allCutsRaw{i}{j} - n;
            comboHeights{j} = (ZID <= 0) * 0 + (ZID > 0) * Z(ZID,3);
        end
        allHeights{i} = comboHeights;
    end

    % 
    allCuts.combos    = allCutsRaw;
    allCuts.leaves    = allLeaves;
    allCuts.heights   = allHeights;
    allCuts.leafCount = allLeafCounts;



    % Nested helper functions


    function cuts = getCutsRestricted(node)

        if node > numel(children) || isempty(children{node})
            if validNodes(node)
                cuts = {{node}};
            else
                cuts = {};
            end
            return;
        end


        cuts = {};
        if validNodes(node)
            cuts{end+1} = {node};
        end


        c = children{node};
        subCombos = {};
        if isscalar(c)
            subCombos = getCutsRestricted(c(1));
        elseif numel(c) == 2
            leftCuts  = getCutsRestricted(c(1));
            rightCuts = getCutsRestricted(c(2));
            for i = 1:length(leftCuts)
                for j = 1:length(rightCuts)
                    subCombos{end+1} = [leftCuts{i}, rightCuts{j}]; %#ok<AGROW>
                end
            end
        end
        cuts = [cuts, subCombos];
    end

    function leaves = getLeaves(node)
        % Recursively collect leaf indices
        if node <= n
            leaves = node;
            return
        end
        if node > numel(childrenPre) || isempty(childrenPre{node})
            leaves = [];
            return
        end
        c = childrenPre{node};
        leaves = [];
        for k = 1:numel(c)
            leaves = [leaves, getLeaves(c(k))]; 
        end
    end
end
