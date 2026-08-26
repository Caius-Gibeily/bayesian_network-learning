function [ss, props, probes, counts, labels, durations] = formatSankey2(clusterArray, propInput,probeInput, epochs)

ss = {};
props = {};
probes = {};
counts = {};
labels = {};

c = 1;
d = 1;

nLayers = size(clusterArray, 2);
durations = [];
for i = 1:nLayers

    clustUnique = unique(clusterArray(:, i));

    % Last layer only produces nodes, no edges
    if i == nLayers
        for j = 1:length(clustUnique)
            clustTraj = clusterArray(clusterArray(:, i) == clustUnique(j), :);

            prop = propInput{i}{j};
            probe = probeInput{i}{j};
            
            probe(prop==0) = [];
            prop(prop==0) = [];
            durations(d) = (epochs(i+1)-epochs(i))/sum(epochs);
            props{d}  = prop;
            probes{d} = probe;
            counts{d} = size(clustTraj, 1);
            labels{d} = sprintf('%d_%d', j, i);
            d = d + 1;
        end
        continue
    end

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

       
        prop = propInput{i}{j};
        probe = probeInput{i}{j};

        probe(prop==0) = [];
        prop(prop==0) = [];
        
        durations(d) = (epochs(i+1)-epochs(i))/sum(epochs);
        props{d}  = prop;
        probes{d} = probe;

        counts{d} = size(clustTraj, 1);
        labels{d} = sprintf('%d_%d', j, i);
        d = d + 1;
    end
end

end
