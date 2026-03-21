function [ROIs2] = localAlign(ROIs1,ROIs2,thresh,maxID)
numROIs1 = unique(ROIs1);
numROIs1 = numROIs1(numROIs1~=0);

numROIs2 = unique(ROIs2);
numROIs2 = numROIs2(numROIs2~=0);

if isempty(numROIs2)
    numROIs2 = numROIs1;
    ROIs2 = ROIs1;
end
overlap = zeros(length(numROIs1),length(numROIs2));

for i = 1:length(numROIs1)
    mask1 = ROIs1 == numROIs1(i);
    for j = 1:length(numROIs2)
        mask2 = ROIs2 == numROIs2(j);
        tt = sum(mask1 == 1 & mask2 == 1,"all");
        %ft = sum(mask1 ~= 1 & mask2 == 1,"all");
        %tf = sum(mask1 == 1 & mask2 ~= 1,"all");
        %ff = sum(mask1 ~= 1 & mask2 ~= 1,"all");
        overlap(i,j) = (tt)/sum(mask2,"all");
    end
end

for i = 1:length(numROIs1)
    [m,idx] = max(overlap(i,:));
    
    if numROIs2(idx) ~= numROIs1(i) && m >= thresh
        occ = numROIs2;
        %slots = 1:max(numROIs2)+1;
        %slot = min(setdiff(slots,occ));
        ROIs2(ROIs2==numROIs1(i)) = maxID+1; %slot
        numROIs2(numROIs2==numROIs1(i)) = maxID+1; %slot
        maxID = max(numROIs2);
        ROIs2(ROIs2==numROIs2(idx)) = numROIs1(i);   
    end


end

