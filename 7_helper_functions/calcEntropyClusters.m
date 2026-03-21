function [H,Hshuf,Hshuf_stats] = calcEntropyClusters(data,clips,reps)

H = cell(size(data));
Hshuf = cell(size(data));
Hshuf_stats = cell(size(data));
for c = 1:length(clips)
    clusterArrayTD = data{c}.clusterArrayTD;
    clusterArrayASD = data{c}.clusterArrayASD;
    
    H{c} = calcEntropy(clusterArrayTD,clusterArrayASD);
    Hshuf{c} = shuffleEntropy(clusterArrayTD,clusterArrayASD,reps);

    Hshuf_stats{c}(:,1) = prctile(Hshuf{c},2.5,2);
    Hshuf_stats{c}(:,2) = median(Hshuf{c},2);
    Hshuf_stats{c}(:,3) = prctile(Hshuf{c},97.5,2);
    
end

function H = calcEntropy(clusterArrayTD,clusterArrayASD)
    H = zeros(size(clusterArrayTD,2),2);
    for i = 1:size(clusterArrayTD,2)

        freqsTD = tabulate(clusterArrayTD(:,i));
        pmf = freqsTD(:,3)/100; pmf(pmf==0) = [];
       
        H(i,1) = -sum(pmf.*log2(pmf))/length(pmf);

        freqsASD = tabulate(clusterArrayASD(:,i));
        pmf = freqsASD(:,3)/100; pmf(pmf==0) = [];

        H(i,2) = -sum(pmf.*log2(pmf))/length(pmf);
    end
end

function Hshuf = shuffleEntropy(clusterArrayTD,clusterArrayASD,reps)
    Hshuf = [];
    nTD = size(clusterArrayTD,1);
    clusterArray = [clusterArrayTD; clusterArrayASD];
    for i = 1:reps  
        idsShuf = randperm(size(clusterArray,1),size(clusterArray,1));
        clusterShuf = clusterArray(idsShuf,:);

        clusterArrayTD_shuf = clusterShuf(1:nTD,:);
        clusterArrayASD_shuf = clusterShuf(nTD+1:end,:);
        
        Hrep = calcEntropy(clusterArrayTD_shuf,clusterArrayASD_shuf);
        Hshuf(:,i) = Hrep(:,1) - Hrep(:,2);
    end

end
        
        
        

end