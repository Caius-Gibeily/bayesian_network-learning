function corrected = correctClusters(clusterArray,indivArray,epochs,tries)

corrected = zeros(size(clusterArray));
for i = 1:size(clusterArray,2)
    clusts = clusterArray(:,i);
    for j = 1:tries
      clusts  = fixClustersBySilhouette(clusts',indivArray(:,epochs(i):epochs(i+1)));
    end
    corrected(:,i) = clusts;
end