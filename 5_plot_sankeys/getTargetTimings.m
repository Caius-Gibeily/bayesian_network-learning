function [on_offset, prop_targetsOb,prop_targetsPerm] = getTargetTimings(targetCell, data, clips, reps)

arguments
    targetCell cell
    data cell
    clips
    reps = 1000
end

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
        display(i)
    end
    on_offset{c}(on_offset{c}(:,1) == 0,:) = [];
end

[on_offset,prop_targetsOb,prop_targetsPerm] = computeProportionsTargets(clips, ...
    data,on_offset,reps);

end