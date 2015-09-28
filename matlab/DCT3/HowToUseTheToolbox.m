%%%%%%%%%%%%%%%% HOW TO USE THE TOOLBOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO STEP COPULA OR COPULA VINE

% STEP ONE (ASSUMES GARCH TYPE MARGINS). CALL:

specM = modelspec(data); %data is the data set. Choose "GARCH model for each series"
% and make the desired choices for the margins. Then call:

[parsM, LogLM, evalM, GradHessM, udata] = fitModel(specM, data, 'fminunc');

% if you do not care for standard errors you can use 'fmincon'. parsM
% contains the estimated parameters, in the form of table 5, in the pdf
% LogLM is the vector of the marginal log likelihoods at the optimum
% evalM is a cell that contains structures
% GradHessM is a cell that contains structures
% udata is the matrix of sample GARCH residuals turned to uniform.
%%
% STEP TWO. 
% Define your copula or copula vine by:
specC = modelspec(udata); 
% Estimate the copula parameters, by:
[parsC, LogLC, evalC, GradHessC] = fitModel(specC, udata, 'fminunc');

% you can use "input defaults" when your are prompt to define the starting
% values. Usually the problem converges to a solution without the use of
% sophisticated initial parameters. Remember that when you type the initial 
% parameters, this should always be a column vector, for example in the DCC
% t copula, the correct way to imput the vector is: [12;.025;.92]
%% 
% STARTING VALUES
% In case you want to estimate the same model in one step it would be wise
% to create a vector of starting values for the one step procedure, as
% follows:
test = [reshape(parsM,[p,1]);parsC]; %p is the dimension
% If you want to fit a copula vine you are advised to use good starting
% values. To obtain such values, one more step is needed, prior to step two
% call:
specS = modelspec(udata); 
% and define the model as "copula vine sequentially"
% Estimate the copula parameters, by:
parsS = fitModel(specS, udata, 'fmincon');
% the array parsS is a good initial guess for the model parameters. Then go
% to step two (line 18)
%%
% ONE STEP COPULA
% Define both margins and copula by:
specM = modelspec(data);
% and selecting "Copula in 1 step", the call:
[pars, LogL, eval, GradHess] = fitModel(spec, udata, 'fminunc');
% for the starting values it would be wise to use the vector "test" created
% in line 33, of the present file
%%

% NOTES:
% 1. Fminunc should be used when asymptotic standard errors are to be
% computed. Fmincon is more robust than fminunc, in the sence that it
% always converges to a solution unlike fminunc

% 2. The estimated parameters are always ordered as in the tables in the
% pdf. For the marginal models, the correct order is:
% [mean parameters;GARCH parameters;marginal distribution parameters]
% For the t DCC copula the order is
% [copula DoF;alphaDCC;bDCC];
% For the SJC copula the first parameter(s) always corresponds to the upper
% tail (equation).




