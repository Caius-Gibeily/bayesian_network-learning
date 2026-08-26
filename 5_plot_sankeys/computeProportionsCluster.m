function [propsTD,propsASD,probesTD,probesASD,probesClustTD,probesClustASD] = computeProportionsCluster(subseqsTD,subseqsASD,on_offset,epoch,clustPatches)

%on_offset = getOnOff(clips,targetCell);
[propsTD,propsASD,probesTD,probesASD,probesClustTD,probesClustASD] = getProps(subseqsTD,subseqsASD,on_offset,epoch,clustPatches);



function on_offset = getOnOff(clips,targetCell)
    
    on_offset = cell(length(clips),1);
    for c = 1:length(clips)
        uVals = cellfun(@(x) unique(x), targetCell{c},'UniformOutput',false);
        mTarget = max(cell2mat(uVals));
        on_offset{c} = zeros(mTarget,3);
        for i = 1:mTarget
            targetPresence = cell2mat(cellfun(@(x) any(x==i), uVals,'UniformOutput',false));
            if any(targetPresence==1) 
                on_offset{c}(i,1) = i;
                [on_offset{c}(i,2),~] = find(targetPresence==1,1,'first');
                [on_offset{c}(i,3),~] = find(targetPresence==1,1,'last');
            end
        end
        on_offset{c}(on_offset{c}(:,1) == 0,:) = [];
    end
end

function [propsTD,propsASD,probesTD,probesASD,probesClustTD,probesClustASD] = getProps(subseqsTD,subseqsASD,on_offset,epoch,clustPatches)
        propsTD = [];
        propsASD = [];

        probes = [];
        probesClust = [];
        on_offsetEp = on_offset;
    
        indiv1 = subseqsTD;
        indiv2 = subseqsASD;
        
        

        startinEpoch = on_offsetEp(:,2) < epoch(2) & on_offsetEp(:,3) > epoch(2);
        on_offsetEp(startinEpoch,3) = epoch(2);

        endinEpoch = on_offsetEp(:,2) < epoch(1) & on_offsetEp(:,3) > epoch(1);
        on_offsetEp(endinEpoch,2) = epoch(1);

        withinEpoch = on_offsetEp(:,3) <= epoch(2) & on_offsetEp(:,2) >= epoch(1);
        on_offsetEp(~withinEpoch,:) = [];
        
        on_offsetEp = [0 epoch(1) epoch(2); on_offsetEp];

        for i = 1:size(on_offsetEp,1)
            probes(i) = on_offsetEp(i,1);
            S1 = indiv1(:,on_offsetEp(i,2):on_offsetEp(i,3));
            S2 = indiv2(:,on_offsetEp(i,2):on_offsetEp(i,3));
         
            propsTD(i) = mean(sum(S1 == on_offsetEp(i,1),1,"omitmissing")./sum(~isnan(S1),1));
            propsASD(i) = mean(sum(S2 == on_offsetEp(i,1),1,"omitmissing")./sum(~isnan(S2),1));
            clustId = table2array(clustPatches(clustPatches.target_id==on_offsetEp(i,1),4));
            if ~isempty(clustId)
                probesClust(i) = clustId;
            else
                probesClust(i) = 0;
            end

        end
        probesASD = probes(propsASD~=0);
        probesClustASD = probesClust(propsASD~=0);
        propsASD(propsASD==0) = [];

        probesTD = probes(propsTD~=0);
        probesClustTD = probesClust(propsTD~=0);
        propsTD(propsTD==0) = [];
        
        

end
        %prop_targets{c}(:,2:3) = prop_targets{c}(:,2:3)./reps;
end

