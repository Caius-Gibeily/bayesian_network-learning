function p_values = computeP(all_logOddsPerm,logOddsOb)

p_values = sum(abs(all_logOddsPerm-0.5) >= abs(logOddsOb-0.5), 2) ./ size(all_logOddsPerm, 2);

end