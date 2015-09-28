% This function plots the contours of likelihood values on the scatter plot
% of a 2 dimensional data.
function [xgrid,ygrid,Z] = biVariateContourPlotsGMMCopula(givenData,gmmObject,logLikelihoodVal,numMeshPoints,x_dim,y_dim)

%INPUT: givenData (MxN, M=number of points, N=Dimension)
%     : plo = binary variable (1 plot contour plot, 0 do not plot)
%OUTPUT: xgrid,ygrid,Z ( Z contains the likelihood values of the points defined by xgrid and ygrid)

%By Ashutosh Tewari (tewaria@utrc.utc.com)
%load general_info;

d = 2;
if nargin < 5
    x_dim = 1;
    y_dim = 2;
end

if x_dim == y_dim
    hist(givenData(:,x_dim),10);
    return;
end

numMeshPoints = min(numMeshPoints,256);

givenData = givenData(:,[x_dim y_dim]);
alpha = gmmObject.alpha;
mu = gmmObject.mu(:,[x_dim y_dim]);
sigma = gmmObject.sigma([x_dim y_dim],[x_dim y_dim],:) + 0.005*repmat(eye(d),[1 1 numel(alpha)]);


gmmObject = gmdistribution(mu,sigma,alpha);

bin_num = 256;
for j = 1:2
   l_limit = min(gmmObject.mu(:,j))-3*(max(gmmObject.Sigma(j,j,:))^0.5);
   u_limit = max(gmmObject.mu(:,j))+3*(max(gmmObject.Sigma(j,j,:))^0.5);
   xmesh_inverse_space{j} = (l_limit:(u_limit-l_limit)/(bin_num-1):u_limit);
end


%if isempty(xmesh)||isempty(pdensity)||isempty(cdensity)
    % Following for loop does the non-parameteric estimation of marginal
    % densities if not provided
    for i = 1:d
        currentVar = givenData(:,i);        
        % finding non-parametric PDF values at prespecified mesh points
        [bandwidth,pdensity{i},xmesh{i}]=kde(currentVar,numMeshPoints);
        pdensity{i}(find(pdensity{i}<0)) = 0;

        % finding non-parametric CDF values
        cdensity{i} = cumsum(pdensity{i});
        cdensity{i} = (cdensity{i}-min(cdensity{i}))/(max(cdensity{i})-min(cdensity{i})); % scaling the cdensity value to be between [0 1]
    end
%end

[xgrid,ygrid] = meshgrid(xmesh{1}(2:end-1),xmesh{2}(2:end-1));

for k = 1:d
    marginalLogLikelihood_grid{k} = log(pdensity{k}(2:end-1)+eps);
    marginalCDFValues_grid{k} = cdensity{k}(2:end-1);
end
[marg1,marg2] = meshgrid(marginalLogLikelihood_grid{1},marginalLogLikelihood_grid{2});

[xg,yg] = meshgrid(marginalCDFValues_grid{1},marginalCDFValues_grid{2});
inputMatrix = [reshape(xg,numel(xg),1) reshape(yg,numel(yg),1)];
clear xg yg;

copulaLogLikelihoodVals = gmmCopulaPDF(inputMatrix,gmmObject,xmesh_inverse_space);
Z = reshape(copulaLogLikelihoodVals,size(marg1,1),size(marg1,2));
Z = Z+marg1+marg2;

Z = exp(Z); % Getting the likelihood value from the log-likelihood.

plot(givenData(:,1),givenData(:,2),'k.','MarkerSize',3);hold
contour(xgrid,ygrid,Z,40);
%title_string = ['GMCM fit (Log-Likelihood = ',num2str(logLikelihoodVal), ')'];
%title(title_string,'FontSize',12,'FontWeight','demi');
axis tight;
