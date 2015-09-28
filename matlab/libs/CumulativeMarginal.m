function pr = CumulativeMarginal( bnData, index )
% marginal probability of variable var in value x
% TODO a kernel density estimation supporting missing values
    data_index = find(~isnan(bnData(:, index)));
    missing_index = setdiff(1:size(bnData, 1), data_index);
    pr = zeros(size(bnData, 1), 1);
    pr(data_index) = ksdensity(bnData(data_index, index), bnData(data_index, index), 'function', 'cdf');
    pr(missing_index) = NaN;
end

