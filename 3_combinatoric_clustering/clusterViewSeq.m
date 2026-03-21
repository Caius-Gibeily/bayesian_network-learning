function clusterArray = clusterViewSeq(indivComb,epochs,k,minCutOff)
clusterArray = zeros(size(indivComb,1),length(epochs)-1);

for i = 1:length(epochs)-1
    window = indivComb(:,epochs(i):epochs(i+1));
    clusterArray(:,i) = getOptimalClusters(window,k,minCutOff);
    display(i)
end