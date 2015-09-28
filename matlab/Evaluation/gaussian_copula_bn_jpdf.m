function pr = gaussian_copula_bn_jpdf( x, y, rho )
%GAUSSIAN_COPULA_BN_JPDF Summary of this function goes here
%   currently only standard gaussian copula is supported, the marginals are
%   estimated by use of ksdensity which is a kernel smooth density
%   estimator. 
%   x is the data point to be evaluated, and y is the dataset used for
%   estimating marginals
    instance = zeros(size(x, 2), 1);
    for i = 1:size(x, 2)
      Fx_i = ksdensity(y(:, i), x(i), 'function', 'cdf');
      instance(i, 1) = norminv(Fx_i);
    end
    mar = find(isinf(instance));
    if ~isempty(mar)
        instance(mar) = [];
        rho(mar, :) = [];
        rho(:, mar) = [];
    end
    mu = zeros(1, length(instance));
    pr = mvnpdf(instance', mu, rho);

end

