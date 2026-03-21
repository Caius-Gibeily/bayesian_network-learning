function allCDFs = deriveCDF(combs,c)


height = zeros(size(combs.heights));
equality = zeros(size(combs.leaves));
nclusts = zeros(size(combs.combos));
for comb = 1:length(combs.combos)
    if isscalar(combs.combos{comb})
        height(comb) = NaN; equality(comb) = NaN; nclusts(comb) = NaN;
    else
        height(comb) = range(cell2mat(combs.heights{comb}));
        nums = cellfun(@(x) length(x), combs.leaves{comb},'UniformOutput',false);
        equality(comb) = range(cell2mat(nums));
    
        nclusts(comb) = abs(length(combs.combos{comb}) - c);
    end
end

hCDF = ecdfIndex(height)'; eCDF = ecdfIndex(equality)'; nCDF = ecdfIndex(nclusts)';
allCDFs = [hCDF eCDF nCDF];


