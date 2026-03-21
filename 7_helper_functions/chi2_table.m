function [chi2_stat, df, p_value, expected] = chi2_table(T)
% CHI2_TABLE  Compute Pearson chi-square statistic for a contingency table.
%
%   [chi2_stat, df, p_value, expected] = chi2_table(T)
%
% INPUT
%   T : r x c contingency table (non-negative counts)
%
% OUTPUT
%   chi2_stat : chi-square statistic
%   df        : degrees of freedom
%   p_value   : p-value (right-tail)
%   expected  : expected counts under independence

    % Totals
    row_totals = sum(T, 2);
    col_totals = sum(T, 1);
    grand_total = sum(row_totals);

    % Expected counts under independence
    expected = (row_totals * col_totals) / grand_total;

    % Chi-square statistic
    chi2_stat = sum((T - expected).^2 ./ expected, 'all');

    % Degrees of freedom
    df = (size(T,1)-1) * (size(T,2)-1);

    % p-value
    p_value = 1 - chi2cdf(chi2_stat, df);
end
