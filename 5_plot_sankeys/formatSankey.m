function [ss, props, probes, counts, labels] = formatSankey(clusterArray, patterns)

ss = {};
props = {};
probes = {};
counts = {};
labels = {};

c = 1;
d = 1;

nLayers = size(clusterArray, 2);

for i = 1:nLayers

    clustUnique = unique(clusterArray(:, i));

    % Last layer only produces nodes, no edges
    if i == nLayers
        for j = 1:length(clustUnique)
            clustTraj = clusterArray(clusterArray(:, i) == clustUnique(j), :);

            prop = tabulate(patterns{i}(j, :));
            prop(prop(:,3)==0,:) = [];

            props{d}  = prop(:,3);
            probes{d} = prop(:,1);
            counts{d} = size(clustTraj, 1);
            labels{d} = sprintf('%d_%d', j, i);
            d = d + 1;
        end
        continue
    end

    % For all but last layer:
    % -------------------------------------------------------
    %  Instead of: j (source cluster) then k (dest)
    %  We now loop: k (dest) then j (source cluster)
    % -------------------------------------------------------
    destStates = unique(clusterArray(:, i+1));

    for k = 1:length(destStates)

        for j = 1:length(clustUnique)

            clustTraj = clusterArray(clusterArray(:, i) == clustUnique(j), :);

            % Compute fraction of transitions from j→dest(k)
            idx = clustTraj(:, i+1) == destStates(k);
            n_trans = sum(idx) / size(clustTraj, 1);

            % Add edge only if nonzero
            if n_trans > 0
                ss{c,1} = char(strcat(string(j), "_{", string(i), "}"));
                ss{c,2} = char(strcat(string(k), "_{", string(i+1), "}"));
                ss{c,3} = n_trans;
                c = c + 1;
            end
        end
    end

    % Node property tables (per cluster j)
    for j = 1:length(clustUnique)
        clustTraj = clusterArray(clusterArray(:, i)==clustUnique(j), :);

        prop = tabulate(patterns{i}(j, :));
        prop(prop(:,3)==0,:) = [];

        props{d}  = prop(:,3);
        probes{d} = prop(:,1);
        counts{d} = size(clustTraj, 1);
        labels{d} = sprintf('%d_%d', j, i);
        d = d + 1;
    end
end

end
