function allCuts = all_cluster_combinations(Z, minLeaves)
    % ALL_CLUSTER_COMBINATIONS Generate all valid cluster combinations after pruning.
    %
    %   allCuts = all_cluster_combinations(Z, minLeaves)
    %
    %   Z         - hierarchical clustering linkage matrix
    %   minLeaves - minimum number of basal leaves a cluster must have
    %
    %   Returns struct:
    %       allCuts.combos : node IDs per combination
    %       allCuts.leaves : leaf IDs per cluster (same shape)

    if nargin < 2
        minLeaves = 5;
    end

    n = size(Z, 1) + 1;               % number of leaves
    nNodes = n + size(Z, 1);
    children = cell(1, nNodes);
    
    % --- Build tree structure ---
    for i = 1:size(Z, 1)
        parent = n + i;
        children{parent} = [Z(i,1), Z(i,2)];
        heights{parent} = [Z(i,3)];
    end
    childrenPre = children;
    % --- Prune nodes with too few leaves ---
    [~, children] = prune(nNodes, children, minLeaves, n);

    allCutsRaw = getCuts(nNodes);

    keyStrings = cell(size(allCutsRaw));
    for i = 1:numel(allCutsRaw)
        combo = allCutsRaw{i};
        if iscell(combo)
            combo = cell2mat(combo);
        end
        keyStrings{i} = sprintf('%d_', sort(combo));
    end
    [~, ia] = unique(keyStrings, 'stable');
    allCutsRaw = allCutsRaw(ia);

    % Compute leaves per cluster
    allLeaves = cell(size(allCutsRaw));
    for i = 1:numel(allCutsRaw)
        combo = allCutsRaw{i};
        if iscell(combo)
            combo = cell2mat(combo);
        end
        leavesPerCluster = cell(size(combo));
        for j = 1:numel(combo)
            leaves = getLeaves(combo(j));
            if isempty(leaves)
                % fallback: if node has no recorded children, treat as leaf
                leaves = combo(j);
            end
            leavesPerCluster{j} = sort(unique(leaves));
        end
        allLeaves{i} = leavesPerCluster;
    end

    % Return structured output

    filt = cellfun(@(x) length(cell2mat(x))==n,allLeaves,'UniformOutput',false);
    allCutsRaw = allCutsRaw(cell2mat(filt)==1);
    allLeaves = allLeaves(cell2mat(filt)==1);

    allCuts.combos = allCutsRaw;
    allCuts.leaves = allLeaves;

    % heights
    allHeights = cell(size(allCutsRaw));
    for comb = 1:length(allCutsRaw)
        for clust = 1:length(allCutsRaw{comb})
            ZID = allCutsRaw{comb}{clust} - n;
            if ZID <= 0
                allHeights{comb}{clust} = 0;
            else
                allHeights{comb}{clust} = Z(ZID,3);
            end
        end
    end
    allCuts.heights = allHeights;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested helper functions 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [nLeaves, children] = prune(node, children, num, n)
        % Recursively prune subtrees with fewer than `num` leaves
        if node <= n
            nLeaves = 1;
            return
        end

        c = children{node};
        if isempty(c)
            nLeaves = 0;
            return
        end

        % Recursively count leaves under children
        [leftLeaves, children] = prune(c(1), children, num, n);
        rightLeaves = 0;
        if numel(c) > 1
            [rightLeaves, children] = prune(c(2), children, num, n);
        end

        nLeaves = leftLeaves + rightLeaves;

        % Remove small subclusters
        if leftLeaves < num
            children{node}(children{node} == c(1)) = [];
        end
        if rightLeaves < num
            children{node}(children{node} == c(2)) = [];
        end
    end


    function cuts = getCuts(node)
        % Recursively list all valid cluster partitions
        if node > numel(children) || isempty(children{node})
            % node is leaf or pruned
            cuts = {{node}};
            return;
        end

        c = children{node};
        cuts = {{node}}; % option 1: treat as whole cluster

        if isscalar(c)
            leftCuts = getCuts(c(1));
            for i = 1:length(leftCuts)
                cuts{end+1} = [leftCuts{i}]; %#ok<AGROW>
            end
        elseif numel(c) == 2
            leftCuts = getCuts(c(1));
            rightCuts = getCuts(c(2));
            for i = 1:length(leftCuts)
                for j = 1:length(rightCuts)
                    cuts{end+1} = [leftCuts{i}, rightCuts{j}]; %#ok<AGROW>
                end
            end
        end
    end


    function leaves = getLeaves(node)
        % Recursively collect leaf indices
        if node <= n
            leaves = node;
            return
        end
        if node > numel(childrenPre) || isempty(childrenPre{node})
            leaves = node * (node <= n); % return itself if leaf-like
            return
        end

        c = childrenPre{node};
        leaves = [];
        for k = 1:numel(c)
            leaves = [leaves, getLeaves(c(k))]; %#ok<AGROW>
        end

        % Keep only valid leaf indices
        leaves = leaves(leaves > 0 & leaves <= n);
    end
end
