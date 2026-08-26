function [ss2,props2,probes2,diffs2,counts2,nodeList,nodeLayer] = filterParent(ss,parent,child,labels,props,probes,diffs,counts)



ssStr = string(ss);

col1 = split(ssStr(:,1), "_");   % N × 2 string array
parentSuffix = col1(:,2);

col2 = split(ssStr(:,2), "_");
childSuffix = col2(:,2);

parent = string(parent);
child  = string(child);

isParent = (parentSuffix == parent);
isChild  = (childSuffix  == child);


ss2 = ss(isParent & isChild,:);
diffs2 = diffs(isParent & isChild,:);


children = ss(isParent & isChild,2);
parents = ss(isParent & isChild,1);
[~, idx] = ismember([parents;children], labels);
nodeList = unique(labels(idx));

parts = split(nodeList', "_");   
suffix = str2double(parts(:,2)); 

[~, order] = sort(suffix);
nodeList = nodeList(order);


for i = 1:length(nodeList)
    mask = string(labels') == nodeList(i);
    props2{i} = props{mask}; 
    probes2{i} = probes{mask};     
    counts2{i} = counts{mask}; 
end

nodeLayer = split(string(nodeList'),"_");
nodeLayer = str2double(nodeLayer(:,2));

end