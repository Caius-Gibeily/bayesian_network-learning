function silsAll = collateSilhouettes(data,clips)

arguments
    data cell
    clips 
end

silsAll = cell(length(clips),2);




for c = 1:length(clips)    
    silsAll{c,1} = collate(data,"TD");
    silsAll{c,2} = collate(data,"ASD");
end


function sils = collate(data,dx)
    sils = [];
    for i = 1:size(data{c}.(strcat("clusterArray",dx)),2)
        epochs = data{c}.epochs;
        dataSub = data{c}.(strcat("indivArray",dx));
        data_clustSub = data{c}.(strcat("clusterArray",dx));
        sils_i = silhouette(dataSub(:,epochs(i):epochs(i+1)),...
        data_clustSub(:,i),'Hamming',"FaceColor","blue");

        sils(i) = mean(sils_i);
    end 
end
end
