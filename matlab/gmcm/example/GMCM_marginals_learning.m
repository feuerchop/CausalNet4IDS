% Fitting the marginal density using a diffusion based non-parameteric density estimation
% approach (refer to the function 'kde' for more details).
function u = GMCM_marginals_learning(givenData,d,numMeshPointVector)

if numel(numMeshPointVector) == 1;
    numMeshPointVector = ones(d,1)*numMeshPointVector;
end

for i = 1:d
    currentVar = givenData(:,i);
    numMeshPoints = numMeshPointVector(i);
    [~,pdensity{i},xmesh{i}]=kde(currentVar,numMeshPoints);

    % finding non-parametric CDF values
    cdensity{i} = cumsum(pdensity{i});
    cdensity{i} = (cdensity{i}-min(cdensity{i}))/(max(cdensity{i})-min(cdensity{i})); % scaling the cdensity value to be between [0 1]
    
    % obtaining the cdf values by linear interpolation
    u(:,i) = interp1q(xmesh{i}',cdensity{i},givenData(:,i));
end
