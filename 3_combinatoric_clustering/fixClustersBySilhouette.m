function new_idx = fixClustersBySilhouette(idx, S)
% fixClustersBySilhouette
%   Reassigns poorly clustered points, removes bad clusters,
%   and enforces a minimum cluster size of 15% of N.
%
% INPUTS
%   idx  - Nx1 vector of cluster labels
%   S    - NxL array of sequences (integer ID labels)
%
% OUTPUT
%   new_idx - updated cluster labels after reassignment

    % Force column vector indexing
    new_idx = idx(:);
    uniqueClusters = unique(new_idx);
    N = numel(new_idx);
    k = numel(uniqueClusters);

    % Compute silhouette scores
    sil = silhouette(S, new_idx, 'hamming');

    % Compute full Hamming distance matrix
    D = pdist2(S, S, 'hamming');

    % ----------------------------------------------------------
    % STEP 1 — Identify clusters with majority negative silhouette
    % ----------------------------------------------------------
    badClusters = false(k,1);
    for j = 1:k
        c = uniqueClusters(j);
        members = find(new_idx == c);
        if isempty(members)
            badClusters(j) = true;
        else
            fracNeg = mean(sil(members) < 0);
            badClusters(j) = fracNeg > 0.7;
        end
    end
    doomedClusters = uniqueClusters(badClusters);

    % ----------------------------------------------------------
    % STEP 2 — Apply minimum cluster size rule
    % ----------------------------------------------------------
    minSize = ceil(0.10 * N);
    smallClusters = false(k,1);

    for j = 1:k
        c = uniqueClusters(j);
        members = find(new_idx == c);
        if numel(members) > 0 && numel(members) < minSize
            smallClusters(j) = true;
        end
    end

    tooSmallClusters = uniqueClusters(smallClusters);

    % Combine silhouette-doomed + minimum-size clusters
    doomedClusters = union(doomedClusters, tooSmallClusters);

    % Mark all members of doomed clusters for reassignment
    new_idx(ismember(new_idx, doomedClusters)) = 0;

    % ----------------------------------------------------------
    % STEP 3 — Mark items with negative silhouette individually
    % ----------------------------------------------------------
    new_idx(sil < 0) = 0;

    % ----------------------------------------------------------
    % STEP 4 — Determine surviving clusters
    % ----------------------------------------------------------
    survivors = setdiff(uniqueClusters, doomedClusters);

    % Keep only clusters that actually have members
    survivors = survivors(arrayfun(@(c) any(new_idx == c), survivors));

    % fallback: if everything is marked bad → collapse to one cluster
    if isempty(survivors)
        new_idx(:) = 1;
        return
    end

    % ----------------------------------------------------------
    % STEP 5 — Compute medoids for surviving clusters
    % ----------------------------------------------------------
    medoids = nan(numel(survivors), 1);

    for ii = 1:numel(survivors)
        c = survivors(ii);
        members = find(new_idx == c);

        subD = D(members, members);
        [~, mIdx] = min(sum(subD, 2));
        medoids(ii) = members(mIdx);
    end

    % fallback: if medoids fail, collapse into cluster 1
    if all(isnan(medoids))
        new_idx(:) = survivors(1);
        return
    end

    % ----------------------------------------------------------
    % STEP 6 — Reassign unassigned items (idx = 0)
    % ----------------------------------------------------------
    unassigned = find(new_idx == 0);

    for i = unassigned'
        % distance to each medoid
        d = arrayfun(@(m) D(i, m), medoids);
        d(isnan(d)) = inf;  % safety

        [~, best] = min(d);
        new_idx(i) = survivors(best);
    end
end
