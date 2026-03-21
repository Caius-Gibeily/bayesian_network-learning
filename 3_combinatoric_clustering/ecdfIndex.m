function ecdf_vals = ecdfIndex(x)

    % Validate input
    if ~isvector(x)
        error('Input X must be a vector.');
    end

    n = numel(x);

    % Compute midranks (handles ties correctly)
    ranks = tiedrank(x);

    % Convert to ECDF values (proportion of samples ≤ x)
    ecdf_vals = ranks / n;
end
