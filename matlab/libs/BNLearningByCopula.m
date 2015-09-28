function [equivDAGs, rho, pcdag, colliders] = BNLearningByCopula( bnData, alpha)
%BNLEARNINGBYCOPULA Summary of this function goes here
%   Detailed explanation goes here
    % number of nodes or variables
    [num_samples, n] = size(bnData);
    dag = zeros(n, n);
    nodes = 1:n;
    marginals = zeros(size(bnData));
    for i = 1:n
        % define nodes fron 1 to n
        marginals(:, i) = CumulativeMarginal(bnData, i);
    end
   % correlation matrix of all variables 
    mCorr = gaussian_copula_fit('Gaussian', marginals);
    % compute covariance ignoring missing values 
%     dataCov = nancov(bnData); 
%     dataCorr = corrcov(dataCov);
    %mInvCorr = inv(mCorr);
    mInvCorr = (mCorr+1e-5*eye(n))\eye(n);
    mScaleInvCorr = ScaleMatrix(mInvCorr);
   
    for i = 1:n
        for j = (i+1):n
            if abs(mScaleInvCorr(j, i)) > alpha
      %      if FisherSignificanceTest(num_samples, n-2, alpha, mScaleInvCorr(j, i))
                dag(i, j) = 1;
                dag(j, i) = 1;
            end
        end
    end
    % guarantee no circle will be introduced by detriangulation
    [dag, colliders] = detriangulation( dag, mCorr, alpha);   
        % Let us maintain a list of undirected edges
    undirectedEdges = [];
    for i = 1:size(dag, 1)
        for j = i+1:size(dag, 2)
            if dag(i, j) == 1 && dag(j, i) == 1
                undirectedEdges = [undirectedEdges; i, j];
            end
        end
    end
    [pcdag, undirectedEdges] = ConstrainPropagation(dag, undirectedEdges);
    equivDAGs = getAllEquivalentClasses(pcdag, {}, undirectedEdges);
    rho = mCorr;
end

