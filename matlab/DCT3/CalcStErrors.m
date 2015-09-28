function [derivatives, RobVCV, VCV, hessian, RobStE]=CalcStErrors(MyFunc, theta, data, grad, hessian, spec, solver)
% Calculate asymptotic standard errors
% Comments:
% The hessian from fminunc is for the unconstrained problem, however my
% problem is constrained. In order to calculate standard errors the hessian
% wrt the constrained parameters is needed. However numerical hessians are
% not always positive definite. This function works as follows:
% First it calculates the numerical hessian wrt the constrained parameters,
% if this hessian is not positive definite the function calculates the "constrained"
% hessian from the "unconstrained" one (that hessian is the one provided by
% the fminunc and is always positive definite) with the chain rule.
% INPUTS: 
% MyFunc:       String with the name of the (likelihood) function
% theta:        vector of actuall (constrained) parameters
% data:         Data matrix for the likelihood
% hessian:      The hessian output from fminunc
% varargin:     additional arguments for the function
% OUTPUT
% derivatives:  Structure array that contains many many things...
% RobVCV:       The robust variance covariance matrix (Godambe information
%               matrix).
% VCV:          Variance covariance matrix (Fisher's information matrix)
% hessian:      The hessian wrt the actuall parameters
% StE:          diagonal elements of VCV. ^(1/2)
% RobStE:       diagonal elements of RobVCV.^(1/2)
% ------------------------------------------------------------------------
% based on hessian_2sided.m function from UCSD_MFE toolbox of Kevin
% Sheppard
% ------------------------------------------------------------------------
% author: Manthos Vogiatzoglou, UoM, 2010
% contact: vogia@yahoo.com
% ------------------------------------------------------------------------
T = size(data,1);
%if strcmp(spec.purpose,'fitCopVine') == 1 && strcmp(spec.family,'SJC')==1
%    n = spec.size;
%    theta = reshape(theta,[n*(n-1),1]);
%end
scores = MyFuncScores(MyFunc,theta,data, spec, 'fmincon');%calculates the gradient
Hnum = hessian_2sided(MyFunc,theta,data,spec,'fmincon'); % calculates the Hessian
% Check if Hnum is positive definite
Eigvals = eig(Hnum);
mineig = min(Eigvals);
if mineig<0 || isreal(Eigvals) == 0
   if isempty(hessian)==1
        error('provide the solver hessian as the fourth argument of CalcStErrors.m');
   else
    display('Numerical hessian is not PSD, standard errors will be calculated with the Hessian from the solver.')
    display(' ');
    if strcmp(solver, 'fminunc')==1
        display('Standard errors are calculated by the chain rule.')
        [dum, firstDer, secDer] = RescaleParameters(theta, 2, spec);
        % calculate the real hessian from the unconstrained one
        Hs = hessian.*(firstDer*firstDer')+diag(grad.*secDer); % the hessian from the chain rule
        VCV = inv(Hs/T); % Fisher's info matrix
        RobVCV = (inv(Hs/T)*cov(scores)*inv(Hs/T))/T; % Godambe info matrix
        RobStE = diag(sqrt(RobVCV));
    elseif strcmp(solver,'fmincon')==1
        display('Standard errors will be calculated with the fmincon hessian. Results might be inaccurate!')
        Hs = hessian;
        VCV = inv(Hs/T); % Fisher's info matrix
        RobVCV = (inv(Hs/T)*cov(scores)*inv(Hs/T))/T; % Godambe info matrix
        RobStE = diag(sqrt(RobVCV));
    end
   end
else
    Hs = Hnum;
    VCV = inv(Hs/T); % Fisher's info matrix
    RobVCV = (inv(Hs/T)*cov(scores)*inv(Hs/T))/T; % Godambe info matrix
    RobStE = diag(sqrt(RobVCV));      
end
derivatives.hessian = Hs;
derivatives.scores = scores;
derivatives.VCV = VCV;
derivatives.RobVCV = RobVCV;
derivatives.RobStE = RobStE;
