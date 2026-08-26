function plotEntropyClips(data,clips)

% Permute group labels and calculate entropy diff
[H,~,Hshuf_stats] = calcEntropyClusters(data,clips,1000);

t = tiledlayout(3,2,'TileSpacing','compact');
for c = 1:length(clips)
    nexttile
    plotEntropy(H{c},Hshuf_stats{c})
end

ylabel(t,"Epoch")
xlabel(t,{"Normalised entropy", ...
    "difference (\DeltaH, TD - ASD)"})

end