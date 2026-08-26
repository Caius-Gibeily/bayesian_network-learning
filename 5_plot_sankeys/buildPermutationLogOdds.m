function [all_log_odds, all_logOddsPerm] = buildPermutationLogOdds(prop_targetsPerm, clips, reps)

n = numel(prop_targetsPerm);

log_odds_cell = cell(n,1);
perm_cell = cell(n,1);

for c = 1:n

    pt = cell2mat(prop_targetsPerm{c});
    id_col = repelem(clips(c), size(pt,1))';

    logOddsPerm = pt(:,2:reps+1) ./ (pt(:,reps+2:end) + pt(:,2:reps+1));

    logOdds_mean = mean(logOddsPerm,2);
    logOdds_low  = prctile(logOddsPerm,2.5,2);
    logOdds_high = prctile(logOddsPerm,97.5,2);

    log_odds_cell{c} = [id_col, logOdds_low, logOdds_mean, logOdds_high];
    perm_cell{c} = logOddsPerm;

end

all_log_odds = vertcat(log_odds_cell{:});
all_logOddsPerm = vertcat(perm_cell{:});

end