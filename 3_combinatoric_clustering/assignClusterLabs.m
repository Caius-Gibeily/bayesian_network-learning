function clustArray = assignClusterLabs(clusterComb)
    n = length(cell2mat(clusterComb));
    clustArray = zeros(n,1);
    for clust = 1:length(clusterComb)
        clustArray(clusterComb{clust}) = clust;
    end

end