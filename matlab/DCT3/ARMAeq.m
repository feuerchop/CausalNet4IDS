function [residuals, mu] = ARMAeq(theta, data, spec)

[T,n] = size(data);
if n > 1
    error('the data should be univariate (one column)')
end
m = spec.mparams;
covs = ones(T,m);
for i = 2:m
    covs(:,i) = [zeros(i-1,1); data(1:end-i+1)];
end

mu = sum(repmat(theta',[T,1]).*covs,2);
residuals = data - mu;

