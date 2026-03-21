function clustCombID = findOptimalClustComb(allCDFs,weightVector)

weighted = allCDFs * weightVector';
[~, clustCombID] = min(weighted);

end