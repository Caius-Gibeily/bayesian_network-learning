function best = findOptimalClusterCombination(combs, targetK, lambda)
    if nargin < 3
        lambda = 1;
    end

    % Defensive check
    if ~isfield(combs, "leafCount") || isempty(combs.leafCount)
        error("combs.leafCount is missing or empty");
    end
    if ~isfield(combs, "combos") || isempty(combs.combos)
        error("combs.combos is missing or empty");
    end

    % Determine valid nodes if not stored
    if isfield(combs, 'validNodes')
        validNodes = combs.validNodes;
    else
        % Assume nodes correspond to 1:numel(leafCount)
        validNodes = 1:numel(combs.leafCount);
    end

    nodeToIndex = containers.Map(validNodes, 1:numel(validNodes));
    bestScore = -inf;
    bestIdx = NaN;
    allScores = nan(1, numel(combs.combos));

    fprintf('Evaluating %d combinations...\n', numel(combs.combos));

    for i = 1:numel(combs.combos)
        combo = combs.combos{i};
        if isempty(combo)
            continue;
        end
        if iscell(combo)
            combo = cell2mat(combo);
        end

        % Map node IDs to indices in leafCount
        idx = [];
        for j = 1:numel(combo)
            if isKey(nodeToIndex, combo(j))
                idx(end+1) = nodeToIndex(combo(j)); %#ok<AGROW>
            end
        end

        if isempty(idx)
            % Print first few skips for debugging
            if i < 10
                fprintf('Skipped combo %d: no valid node matches.\n', i);
            end
            continue;
        end

        leafCounts = combs.leafCount(idx);
        Kc = numel(combo);
        varL = var(leafCounts);
        meanL = mean(leafCounts);

        score = -abs(Kc - targetK) - lambda * varL / (meanL + eps);
        allScores(i) = score;

        if score > bestScore
            bestScore = score;
            bestIdx = i;
        end
    end

    % Find next-best if none met target
    if isnan(bestIdx) || isinf(bestScore)
        fprintf('No exact match found, looking for next best...\n');
        validScores = find(~isnan(allScores) & ~isinf(allScores));
        if isempty(validScores)
            warning('No valid cluster combination found.');
            best = struct('idx', [], 'combo', [], 'leaves', [], ...
                          'leafCount', [], 'score', -inf);
            return;
        else
            [~, relIdx] = max(allScores(validScores));
            bestIdx = validScores(relIdx);
            bestScore = allScores(bestIdx);
        end
    end

    best.idx = bestIdx;
    best.combo = combs.combos{bestIdx};
    best.leaves = combs.leaves{bestIdx};

    if iscell(best.combo)
        comboNumeric = cell2mat(best.combo);
    else
        comboNumeric = best.combo;
    end

    idx = [];
    for j = 1:numel(comboNumeric)
        if isKey(nodeToIndex, comboNumeric(j))
            idx(end+1) = nodeToIndex(comboNumeric(j)); %#ok<AGROW>
        end
    end

    best.leafCount = combs.leafCount(idx);
    best.score = bestScore;
    best.Kc = numel(comboNumeric);
    best.varL = var(best.leafCount);
    best.meanL = mean(best.leafCount);

    fprintf('✅ Best combination: %d clusters, mean=%.2f, var=%.2f, score=%.3f\n', ...
        best.Kc, best.meanL, best.varL, best.score);
end
