function pattern = findClustStrategy(data)

clips = cell2mat(cellfun(@(x) x.clip, data, 'UniformOutput', false));
clusterArrayTD = cellfun(@(x) x.clusterArrayTD, data, 'UniformOutput', false);
clusterArrayASD = cellfun(@(x) x.clusterArrayASD, data, 'UniformOutput', false);



pattern = cell(length(clips),1);
for c = 1:length(clips)
    seqsTD = data{c}.recoded.TD;
    seqsASD = data{c}.recoded.ASD;
    % 
    %seqsTD = data{c}.indivArrayTD;
    %seqsASD = data{c}.indivArrayASD;

    epochs = data{c}.epochs;
    
    pattern{c}.TD = getPatterns(epochs,seqsTD,clusterArrayTD{c});
    pattern{c}.ASD = getPatterns(epochs,seqsASD,clusterArrayASD{c});
    

end

end

function props = getProps(epochs,seqsTD,seqsASD,clusterArray, on_offset)
    
    for n = 1:length(epochs)-1
        
        [on_offset,prop_targetsOb,prop_targetsPerm] = computeProportionsCluster(seqsTD,seqsASD,on_offset,reps)


function pattern = getPatterns(epochs,seqs,clusterArray)
    pattern = cell(1,length(epochs)-1);
    for i = 1:length(epochs)-1
        subseqs = seqs(:,epochs(i):epochs(i+1));
        clusts = unique(clusterArray(:,i));
        pattern{i} = zeros(length(clusts),size(subseqs,2));

        for j = 1:length(clusts)

            parts = subseqs(clusterArray(:,i)==clusts(j),:);

            [~,pattern{i}(j,:),~,~] = kmedoids(parts,1);%mode(parts,1);

        end

    end
end





        