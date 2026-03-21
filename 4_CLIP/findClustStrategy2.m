function [props,probes,probesClust] = findClustStrategy2(data,clips,on_offset,clust_identity)

probes = cell(length(data),1);
props = cell(length(data),1);
probesClust = cell(length(data),1);
for c = 1:length(clips)
    %seqsTD = data{c}.recoded.TD;
    %seqsASD = data{c}.recoded.ASD;
    clusterArrayTD = data{c}.clusterArrayTD;
    clusterArrayASD = data{c}.clusterArrayASD;
    % 
    seqsTD = data{c}.indivArrayTD;
    seqsASD = data{c}.indivArrayASD;

    epochs = data{c}.epochs;
    
    [props{c},probes{c},probesClust{c}] = getProps(epochs,seqsTD,seqsASD, ...
        clusterArrayTD,clusterArrayASD, on_offset{c},clust_identity(clust_identity.clip_id==clips(c),:));


end
end

function [props,probes,probesClust] = getProps(epochs,seqsTD,seqsASD, ...
    clusterArrayTD,clusterArrayASD, on_offset,clustPatches)

    props = struct;
    probes = struct;
    probes.TD = cell(length(epochs)-1,1);
    probes.ASD = cell(length(epochs)-1,1);

    probesClust.TD = cell(length(epochs)-1,1);
    probesClust.ASD = cell(length(epochs)-1,1);

    props.TD = cell(length(epochs)-1,1);
    props.ASD = cell(length(epochs)-1,1);

    for i = 1:length(epochs)-1

        
        %subseqsTD = seqsTD(:,epochs(i):epochs(i+1));
        %subseqsASD = seqsASD(:,epochs(i):epochs(i+1));
        
        clusts = unique([clusterArrayTD(:,i);clusterArrayASD(:,i)]);
      
        

        for j = 1:length(clusts)
            
            partsTD = seqsTD(clusterArrayTD(:,i)==clusts(j),:);
            partsASD = seqsASD(clusterArrayASD(:,i)==clusts(j),:);
            [props.TD{i}{j},props.ASD{i}{j},probes.TD{i}{j},probes.ASD{i}{j},probesClust.TD{i}{j},probesClust.ASD{i}{j}] = computeProportionsCluster(partsTD, ...
                partsASD,on_offset,epochs(i:i+1),clustPatches);
        end

    end
end









        