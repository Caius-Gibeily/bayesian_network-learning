function recoded = interpolateCluster(embeddings, data)


clips = cell2mat(cellfun(@(x) x.clip, data, 'UniformOutput', false));
indivArrayTD = cellfun(@(x) x.indivArrayTD, data, 'UniformOutput', false);
indivArrayASD = cellfun(@(x) x.indivArrayASD, data, 'UniformOutput', false);

recoded = cell(length(clips),1);

for c = 1:length(clips)
    
    embed_filt = embeddings(embeddings.clip_id==clips(c),:);

    recoded{c}.TD = interpolate(indivArrayTD{c},embed_filt);
    recoded{c}.ASD = interpolate(indivArrayASD{c},embed_filt);
    display(c)
end
end

function recoded = interpolate(indivArray,embed_filt)
    recoded = zeros(size(indivArray));
    for i = 1:size(indivArray,1)
        
        for j = 1:size(indivArray,2)
            target = indivArray(i,j);
            target_filt = embed_filt(embed_filt.target_id==target,:);
            if isempty(target_filt)
                recoded(i,j) = target;
            else

                [~, idx] = min(abs(target_filt.frame - j));
                %[~, probeMax] = max(table2array(embed_filt(idx,7:end-1)));
                recoded(i,j) = embed_filt.clusts(idx); %probeMax;
            end
        end
    end
end
                
