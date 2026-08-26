function [diff1,diff2] = subtractMatchingPairs(A, B)

% A, B are cell arrays:
% {src, dst, value}
% src/dst formatted like X_{Y}

%% ---- BUILD MAPS ----
Amap = containers.Map;
Bmap = containers.Map;

for i = 1:size(A,1)
    key = sprintf('%s|%s', A{i,1}, A{i,2});
    Amap(key) = A{i,3};
end

for i = 1:size(B,1)
    key = sprintf('%s|%s', B{i,1}, B{i,2});
    Bmap(key) = B{i,3};
end

%% ---- BUILD OUTPUT USING UNION OF KEYS ----
keysA = Amap.keys;
keysB = Bmap.keys;
allKeys = unique([keysA keysB]);

tmp1 = {};
tmp2 = {};
c = 1;
d = 1;

for i = 1:length(allKeys)
    key = allKeys{i};
    parts = split(key,'|');


    if isKey(Amap,key) && isKey(Bmap,key)
        tmp1{c,1} = parts{1};   % src
        tmp1{c,2} = parts{2};   % dst


        tmp2{d,1} = parts{1};   % src
        tmp2{d,2} = parts{2};   % dst
        

        tmp1{c,3} = Amap(key) - Bmap(key);   % matched
        tmp2{d,3} = Amap(key) - Bmap(key);   % matched
        c = c + 1;
        d = d + 1;
        
    elseif isKey(Amap,key)
        tmp1{c,1} = parts{1};   % src
        tmp1{c,2} = parts{2};   % dst

        tmp1{c,3} = Amap(key);   
        c = c + 1;
    else
        tmp2{d,1} = parts{1};   % src
        tmp2{d,2} = parts{2};   % dst
        
        tmp2{d,3} = -Bmap(key);  
        d = d+1;
    end

end

diff1 = sortSrcDst(tmp1);
diff2 = sortSrcDst(tmp2);

function C = sortSrcDst(tmp)
    extractNums = @(s) sscanf(s,'%d_{%d}');
    
    srcNums = cellfun(@(s) extractNums(s)', tmp(:,1), 'UniformOutput', false);
    dstNums = cellfun(@(s) extractNums(s)', tmp(:,2), 'UniformOutput', false);
    
    srcMat = cell2mat(srcNums);   % Nx2 -> [X Y]
    dstMat = cell2mat(dstNums);   % Nx2 -> [X Y]
    
    sortMat = [srcMat dstMat];
    [~, idx] = sortrows(sortMat,[4 3]);
    
    C = tmp(idx,:);
end

end
