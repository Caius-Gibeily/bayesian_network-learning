function plotTargetLogOdds(prop_targets_all, all_log_odds, all_logOddsPerm, embeddings_probe)

figure; hold on

% Observed log-odds
logOddsOb = prop_targets_all(:,3) ./ ...
           (prop_targets_all(:,4) + prop_targets_all(:,3));

[x_sorted, order] = sort(logOddsOb,"descend");
x = 1:numel(x_sorted);

% Colormap
cmap = parula(5);

% Cluster identities
clust_identity = groupsummary(embeddings_probe, ...
                              ["clip_id","target_id"], ...
                              'mode','clusts');

clusters = clust_identity.mode_clusts(order);

% Observed curve
plot(x, x_sorted,'LineWidth',3,'Color','black')

% Permutation mean
plot(x, all_log_odds(order,3),'o-', ...
    'MarkerFaceColor','black', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize',5)

% Compute significance
p_values = computeP(all_logOddsPerm,logOddsOb);
p_adj = mafdr(p_values);

sig = p_adj < 0.05;

% Significance markers
scatter(x(sig),0.85*ones(sum(sig),1), ...
        100,clusters(sig),"|")

% Colored observed markers
scatter(x(sig),x_sorted(sig), ...
        40,clusters(sig),'filled')

% Formatting
colormap(cmap)
clim([1 5])

xlabel("All movie targets (sorted)")
ylabel("Percentage difference, % (TD/(ASD + TD))")
xlim([0 length(x)])

saveas(gcf,"output/target_props.svg")

end