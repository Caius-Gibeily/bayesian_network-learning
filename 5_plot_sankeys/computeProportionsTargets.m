function [on_offset,prop_targetsOb,prop_targetsPerm] = computeProportionsTargets(clips,data,on_offset,reps)
if nargin == 2
    reps = 1000;
end
%on_offset = getOnOff(clips,targetCell);
prop_targetsOb = getProps(on_offset,data,"observed",1);
prop_targetsPerm = getProps(on_offset,data,"permuted",reps);



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

function prop_targets = getProps(on_offset,data,mode,reps)
    prop_targets = {};
    for c = 1:length(clips)
        prop_targets{c} = cell(size(on_offset{c},1),3);
    
        
        for rep = 1:reps
            if mode == "permuted"
                comb = [data{c}.indivArrayTD; data{c}.indivArrayASD];
                labels = [repelem(0,size(data{c}.indivArrayTD,1)), repelem(1,size(data{c}.indivArrayTD,1))];
                labels = labels(randperm(length(labels),length(labels)));
                
                indiv1 = comb(labels==0,:);
                indiv2 = comb(labels==1,:);
            
        
            elseif mode == "observed"
                indiv1 = data{c}.indivArrayTD;
                indiv2 = data{c}.indivArrayASD;
                
            end
            for i = 1:size(on_offset{c},1)
                prop_targets{c}{i,1} = on_offset{c}(i,1);
                S1 = indiv1(:,on_offset{c}(i,2):on_offset{c}(i,3));
                S2 = indiv2(:,on_offset{c}(i,2):on_offset{c}(i,3));
             
                prop_targets{c}{i,2}(rep) = mean(sum(S1 == on_offset{c}(i,1),1)./sum(S1 > 0,1));
                prop_targets{c}{i,3}(rep) = mean(sum(S2 == on_offset{c}(i,1),1)./sum(S2 > 0,1));
        
            end
            display(rep)
        end
        %prop_targets{c}(:,2:3) = prop_targets{c}(:,2:3)./reps;
    end

end

end
