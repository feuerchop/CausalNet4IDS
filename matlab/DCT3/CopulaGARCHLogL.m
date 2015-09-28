function [LogL, LL, varargout] = CopulaGARCHLogL(theta,data,spec,solver)

if nargin == 3
    solver = 'fmincon';
end
if strcmp(solver,'fminunc')==1
    % the input vector theta should be unconstrained
    % but the log likelihood accepts only constrained parameters
    %theta = RescaleParameters(theta, 1, spec);
    % thus the initial theta is transformed to the corresponding
    % constrained theta
end
T = size(data,1);
n = spec.size;
m = spec.vecsize;
mtheta = theta(1:n*m); mtheta = reshape(mtheta,[m,n]);
ctheta = theta(n*m+1:end);
mLogL = zeros(1,n); mLL = zeros(T,n); ht = mLL; stdRes = mLL; udata = mLL;
for i = 1:n
    
    [mLogL(i), mLL(:,i), ht(:,i), stdRes(:,i), udata(:,i)]=GARCHLogL(mtheta(:,i),data(:,i),spec,solver);
end

[cLogL,cLL,tvpar]=CopulaLogL(ctheta,udata,spec, solver);
varargout{1} = tvpar;
LogL = sum(mLogL) + cLogL;

LL = cLL + sum(mLL,2);

