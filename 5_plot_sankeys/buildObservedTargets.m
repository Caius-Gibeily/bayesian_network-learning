function prop_targets_all = buildObservedTargets(prop_targetsOb, clips)

n = numel(prop_targetsOb);
out = cell(n,1);

for c = 1:n
    pt = cell2mat(prop_targetsOb{c});
    id_col = repelem(clips(c), size(pt,1))';
    out{c} = [id_col, pt];
end

prop_targets_all = vertcat(out{:});

end