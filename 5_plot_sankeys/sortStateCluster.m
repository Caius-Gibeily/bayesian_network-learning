function C = sortStateCluster(tmp)

extractNums = @(s) sscanf(s,'%d_{%d}');

srcNums = cellfun(@(s) extractNums(s)', tmp(:,1), 'UniformOutput', false);

srcMat = cell2mat(srcNums);   % Nx2 -> [X Y]


[~, idx] = sortrows(srcMat,[2,1]);
C = tmp(idx,:);
end