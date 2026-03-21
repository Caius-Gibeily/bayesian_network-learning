function clusterLabs = getOptimalClusters(seqArray,k,minCutOff)

distMatrix = squareform(pdist(seqArray, 'hamming')); 
Z = linkage(distMatrix, 'ward');

combs = getclustCombinations(Z,minCutOff);

for i = 1:length(combs.combos)
    if ~isbetween(length(combs.combos{i}),[k-1,k+1])
        combs.combos{i} = [];
        combs.leaves{i} = [];
        combs.heights{i} = [];
    end
end
combs.combos(cellfun(@isempty,combs.combos))=[];
combs.leaves(cellfun(@isempty,combs.leaves))=[];
combs.heights(cellfun(@isempty,combs.heights))=[];

if isempty(combs.combos)
    clusterLabs = ones(size(seqArray,1),1);
else
    
    allCDFs = deriveCDF(combs,k); %ID = 1:size(allCDFs,1);
    % allCDFs = [ID' allCDFs]
    
    % Sort by clust size then equality
    weightVector = [0.01 0.001];
    clustCombID = findOptimalClustComb(allCDFs(:,2),weightVector(2));
    clusterLabs = assignClusterLabs(combs.leaves{clustCombID});
end


end