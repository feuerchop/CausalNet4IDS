% This function plots the contours of likelihood values on the scatter plot
% of a 2 dimensional data.
function [xgrid,ygrid,likelihood_mat_mesh] = biVariateContourPlotsGMM(givenData,gmmObject,numMeshPoints,x_dim,y_dim)

%INPUT: givenData (MxN, M=number of points, N=Dimension)
%     : plo = binary variable (1 plot contour plot, 0 do not plot)
%OUTPUT: xgrid,ygrid,Z ( Z contains the likelihood values of the points defined by xgrid and ygrid)

%By Ashutosh Tewari (tewaria@utrc.utc.com)

d = 2;

if nargin < 4
    x_dim = 1;
    y_dim = 2;
end

if x_dim == y_dim
    hist(givenData(:,x_dim),10);
    return;
end

numMeshPoints = min(numMeshPoints,256);

givenData = givenData(:,[x_dim y_dim]);
PComponents = gmmObject.PComponents;
mu = gmmObject.mu(:,[x_dim y_dim]);
Sigma = gmmObject.Sigma([x_dim y_dim],[x_dim y_dim],:) + 0.005*repmat(eye(d),[1 1 numel(PComponents)]);

gmmObject = gmdistribution(mu,Sigma,PComponents);

logLikelihoodVal = sum(log(pdf(gmmObject,givenData)));

% Following for loop does the non-parameteric estimation of marginal
% densities.
for i = 1:d
    currentVar = givenData(:,i);
%     binWidth = 2*iqr(currentVar)/(numel(currentVar))^(1/3); % bin width based on Freedman�Diaconis rule
%     if binWidth == 0 ; binWidth = range(currentVar)/50; end
%     numMeshPoints = round(range(currentVar)/binWidth);         
    % finding non-parametric PDF values at prespecified mesh points
    [~,~,xmesh{i}]=kde(currentVar,numMeshPoints);

end

[xgrid,ygrid] = meshgrid(xmesh{1},xmesh{2});
inputMatrix = [reshape(xgrid,numel(xgrid),1) reshape(ygrid,numel(ygrid),1)];
likelihood_vector_mesh = pdf(gmmObject,inputMatrix);
likelihood_mat_mesh = reshape(likelihood_vector_mesh,size(xgrid,1),size(ygrid,2));

% figure;
% subplot(1,2,1);
plot(givenData(:,1),givenData(:,2),'kx','MarkerSize',3);hold
contour(xgrid,ygrid,likelihood_mat_mesh,40);

title_string = ['GMM Fit (Log-Likelihood = ',num2str(logLikelihoodVal), ')'];    
title(title_string,'FontSize',12,'FontWeight','demi');

% idx_gmm = cluster(gmmObject,givenData);

% subplot(1,2,2)
% cluster_scatter_plot(givenData,size(gmmObject.Sigma,3),idx_gmm);
