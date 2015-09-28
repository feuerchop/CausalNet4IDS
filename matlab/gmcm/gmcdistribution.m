% Author: Ashutosh Tewari (tewaria@utrc.utc.com)
% Affiliation: Decision Support & Machine Intelligence
%              United Technologies Research Center
%              East Hartford, CT: 06118
% Code based on the paper 'Parametric Characterization of MultiModal Densities with non-Gaussian
% Nodes', A. Tewari, M.J. Giering, A. Raghunathan, OEDM workshop, ICDM 2011 (paper included in the package) 

% DEFINING GMC (Gaussian Mixture Copula) CLASS
classdef gmcdistribution
    
    properties (SetAccess = private)
        mu     % k x d mean vectors
        sigma  % d x d x k  covariance matrices
        alpha  % k x 1  mixing proportions
    end
    
    methods
        % CONSTRUCTOR METHOD
        % Constructs a gmcdistribution object for a given parameter set.
        function obj = gmcdistribution(mu,sigma,alpha)  
            obj.mu = mu;
            obj.sigma = sigma;
            obj.alpha = alpha;
        end
        
        % METHOD TO FIT A GMC DISTRIBUTION TO DATA
        % Inputs: u  :N x d matrix of CDF values obtained after fitting marginal densities.
        %       : K  : number of clusters
        %       : d  : number of dimensions
        %       : N  : number of data samples.
        % varargin:  'Start' :  'rand_init' OR 'EM_init'  (default'EM_init');
        %            'replicates' : Integer  (default 1)
        %            'iteration'  : Integer  (default 100)
        %            'algorithm'  : 'active-set' OR 'interior-point' (default 'active-set')
        % Example:
        % obj.fit(u,K,d,N,'Start','rand_init','replicates',20,'algorithm','interior-point')
        
        function obj = fit(obj,u,K,d,N,varargin)
           
            % parsing the input argument
            p = inputParser;
            p.addParamValue('Start','EM_init');
            p.addParamValue('replicates','1');
            p.addParamValue('iteration','100');
            p.addParamValue('algorithm','active-set');
            p.parse(varargin{:});

            % Choosing the initialization method. the 'rand_init' option
            % randomly assigns the parameters. 'EM_init' option usually generate
            % better inital guess (refer to the paper mentioned above).
            if strcmp(p.Results.Start,'rand_init')
                [x_init,bounds,A,b] = initializeParameters(d,K);
            elseif strcmp(p.Results.Start,'EM_init')
                [x_init,bounds,A,b] = initializeParameters(d,K);
                [x_init,~] = gmcm_EM(u,K,d,x_init);
            end
            
            % Performing the nonlinear optimization
            [x_final,~] = gmcmOptimization(u,K,d,N,bounds,A,b,x_init,str2double(p.Results.iteration),p.Results.algorithm,str2double(p.Results.replicates));
            
            % converting the paramters in vectorized form into arrays
            [mu,sigma,alpha] = vector2GMMParameters(x_final,d,K);
            
            % updating the propoerties of the gmcdistribution object
            obj.mu = mu;
            obj.sigma = sigma;
            obj.alpha = alpha;           
            
        end
        
        
        %METHOD TO SAMPLE FROM THE GMC DISTRIBUTION 
        % Inputs: obj = GMC object
        %           N = number of samples required;
        % Output: samples = N x d matrix where d is the data dimension
        function gmc_samples = random(obj,N)
            
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.sigma,obj.alpha);
            % Sampling from the gmm object
            gmm_samples = random(gmmObject,N);
            
            % Obtaining the marginal distribution of the gmm.
            gmm_marginals = obtainMarginalsOfGMM(obj.mu,obj.sigma,obj.alpha,K,d);
            
            % samples from the gmc are nothing but the marginal cdf values
            % of the gmm samples.
            gmc_samples = nan(size(gmm_samples));
            for i=1:d
                gmc_samples(:,i) = cdf(gmm_marginals{i},gmm_samples(:,i));
            end
            
        end
        
        % METHOD TO CLUSTER THE DATA GIVEN A GMC OBJECT
        % Inputs: obj = GMC object
        %         u   = N x d data to be clustered
        % Output: idx = N x 1 vector of cluster indices.
        function idx = cluster(obj,u)
            
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            N = size(u,1);
            
            % Obtaining the inverse values with respect to the GMM marginal
            % distributions
            inverseVals = computeInverseVals_vectorized(obj.mu,obj.sigma,obj.alpha,u,d,K,N,0);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.sigma,obj.alpha);
            
            % Cluserting the inverse values using the gmm object.
            idx = cluster(gmmObject,inverseVals);
        end
        
        % METHOD TO COMPUTE THE PDF VALUES W.R.T THE GMC DISCTRIBUTION
        % OBJECT
        %Inputs:  obj = GMC object
        %           u = N x d data to be clustered
        %Output: pdf_vals = N x 1 vector of the pdf values.
        function pdf_vals = pdf(obj,u)
        
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            N = size(u,1);
            
            % Obtaining the inverse values with respect to the GMM marginal
            % distributions
            inverseVals = computeInverseVals_vectorized(obj.mu,obj.sigma,obj.alpha,u,d,K,N,0);
                       
            % Obtaining the log-likelihood of the numerator
            small_mat = 1e-323*ones(N,1);
            first_part = zeros(size(inverseVals,1),K);
            inverseVals_hat = zeros(N,d,K);
            for k = 1:K
                V_mat = chol(obj.sigma(:,:,k))';
                inverseVals_hat(:,:,k) = inverseVals - repmat(obj.mu(k,:),N,1); % Getting the mean adjusted inverse vals
                temp_mat = inverseVals_hat(:,:,k)*(inv(V_mat))';
                first_part(:,k) = obj.alpha(k)*(1/(2*pi)^(d/2))*(1/prod(diag(V_mat)))*exp(-0.5*sum(temp_mat.*temp_mat,2));
                clear temp_mat;
            end
            first_part_ll = log(sum(first_part,2) + small_mat);  % A small positive number is added to avoid log(0)
            clear inverseVals;

            % Getting the log-likelihood of the denominator
            second_part = zeros(N,K);
            for j = 1:d
                temp_vector = zeros(N,K);
                for k = 1:K
                    temp_vector(:,k) =  obj.alpha(k)*(1/sqrt(2*pi*obj.sigma(j,j,k)))...
                        *exp(-0.5*(1/obj.sigma(j,j,k))*(inverseVals_hat(:,j,k).^2));
                end
                second_part(:,j) = log(sum(temp_vector,2)+ small_mat);
            end
            second_part_ll = sum(second_part,2);
            clear inverseVals_hat;
            
            log_likelihood = first_part_ll - second_part_ll;
            pdf_vals = exp(log_likelihood);
        
        end
        
        % METHOD TO COMPUTE THE CDF VALUES W.R.T THE GMC DISCTRIBUTION
        % OBJECT
        %Inputs:  obj = GMC object
        %           u = N x d data to be clustered
        %Output: cdf_vals = N x 1 vector of the pdf values.
        function cdf_vals = cdf(obj,u)
        
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            N = size(u,1);
            
            % Obtaining the inverse values with respect to the GMM marginal
            % distributions
            inverseVals = computeInverseVals_vectorized(obj.mu,obj.sigma,obj.alpha,u,d,K,N,0);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.sigma,obj.alpha);
            
            cdf_vals = cdf(gmmObject,inverseVals);
            
        end
    
        
    end
end