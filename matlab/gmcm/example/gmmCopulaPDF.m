% This function model the joint copula density based on a Gaussian Mixture Density

%INPUT: inputVector (uniformly distributed in [0 1] )
%       gmmObject = the MATLAB object of the trained GMM
%       xmesh = Cell array containing x-grid for each dimension
%OUTPUT: logLikelihoodVal (log-likelihood value of the inputVector based on gmm-copula)

%By Ashutosh Tewari (tewaria@utrc.utc.com)

function logLikelihoodVal = gmmCopulaPDF(inputVector,gmmObject,xmesh)

inputVector(inputVector == 0) = inputVector(inputVector == 0) + eps;
inputVector(inputVector == 1) = inputVector(inputVector == 1) - eps;

numModes = gmmObject.NComponents;
numDimensions =  gmmObject.NDimensions;

% Obtaining the Marginals of the GMM.
for i = 1:numDimensions
    for j = 1:numModes
        mu(j,1)  = gmmObject.mu(j,i);
        Sigma(:,:,j) = gmmObject.Sigma(i,i,j);
        p(j) = gmmObject.PComponents(j);
    end
    
    marginalGMM{i} = gmdistribution(mu,Sigma,p);
    
end

% Getting the inverse values of univariate gmm CDFs empirically
for i = 1:numDimensions
    cdfVals{i} = cdf(marginalGMM{i},xmesh{i}');
    inverseVals(:,i) = inverseCDF(inputVector(:,i),xmesh{i},cdfVals{i}');
end

copulaDensityNumerator = pdf(gmmObject,inverseVals)+  1e-323; %Adding a small number if the density value comes out to be zero

copulaDensityDenominator = 0;
for i = 1:numDimensions
    copulaDensityDenominator = copulaDensityDenominator + log(pdf(marginalGMM{i},inverseVals(:,i))+  1e-323);
end

logLikelihoodVal = log(copulaDensityNumerator) - copulaDensityDenominator;



% This function numerically calculated the inverse of a CDF using linear
% interpolation.
function x = inverseCDF(u,xVals,cdfVals)

%adding small white noise to cdfVals to make them unique for interpolation
while(1)
%    cdfVals = cdfVals + sort(wgn(1,numel(cdfVals),rand(1)*eps,rand(1)*eps));
     cdfVals = cdfVals + sort(abs(random('norm',0,10^-8,numel(cdfVals),1)))';

    if numel(cdfVals) == numel(unique(cdfVals))
        break;
    end
end
%cdfVals = cdfVals + sort(10^-10*randn(1,numel(cdfVals)));

x = interp1(cdfVals,xVals,u,'cubic');   
    