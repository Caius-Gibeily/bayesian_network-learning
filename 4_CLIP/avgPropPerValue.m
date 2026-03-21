function [avgProp,avgProbe] = avgPropPerValue(clips,props,probesClust)

avgProp = cell(length(clips),1);
avgProbe = cell(length(clips),1);

for c = 1:length(clips)
    [avgProp{c}.TD,avgProbe{c}.TD] = computeAvg(props{c},probesClust{c},"TD");
    [avgProp{c}.ASD,avgProbe{c}.ASD] = computeAvg(props{c},probesClust{c},"ASD");
end


function [avgProp,avgProbe] = computeAvg(props,probesClust,dx)
        
    probeCells = probesClust.(dx);   
    propCells  = props.(dx);        
    
    nI = numel(probeCells);           
    
    outC = cell(size(probeCells));    
    probesC = cell(size(probeCells));
    for i = 1:nI
        nJ = numel(probeCells{i});    
        
        for j = 1:nJ
            
            vals   = probeCells{i}{j};    
            td = propCells{i}{j};    
            
            avgVec = zeros(1,6);
            
            for v = 0:5
                mask = (vals == v);
                if any(mask)
                    avgVec(v+1) = mean(td(mask));
                else
                    avgVec(v+1) = 0;   
                end
            end
            
            outC{i}{j} = avgVec;
            outC{i}{j} = outC{i}{j}./sum(outC{i}{j}); 
            probesC{i}{j} = 0:5;
        end
    end
    avgProbe = probesC;
    avgProp = outC;              
end

end