function [datIntraClust,datInterClust] = testIntraclustSim(clusters,distMatrix,nperms)

datIntraClust = cell(length(unique(clusters)),1);
datInterClust = cell(length(unique(clusters)),1);

[~,clustIds] = groupcounts(clusters);

for i = 1:length(clustIds)
    datIntraClust{i} = [];
    datInterClust{i} = [];

    [ids,~] = find(clusters==clustIds(i));
    [nonids,~] = find(clusters~=clustIds(i));
    for j = 1:nperms
        testIn = randsample(ids,2);
        testOut = randsample(nonids,1);
        datIntraClust{i} = [datIntraClust{i}, distMatrix(testIn(1),testIn(2))];
        outer = randsample(1:2,1);
        datInterClust{i} = [datInterClust{i}, distMatrix(testIn(outer),testOut(1))];
        %datInterClust{i} = [datInterClust{i}, distMatrix(testIn(2),testOut(1))];
    end
    

        
end