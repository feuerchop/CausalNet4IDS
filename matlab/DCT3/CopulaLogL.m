function [LogL,LL,varargout]=CopulaLogL(theta,data,spec, solver)
% ---- Log Likelihood functions of the supported copulas -----
% INPUTS:
% theta:        vector of parameters
% data:         matrix with U(0,1) margins
% CopulaSpec:   structured array that contains the various input arguments
%               that define the model. To obtain run the function
%               setCopulaLLinputs.m
% OUTPUTS:       
% LogL:         The negative log - likelihood of the corresponding copula
% LL:           array of negative log likelihoods at each data point
% Rt:           The evolution of the copula parameter
% ------------------------------------------------------------------------
% author: Manthos Vogiatzoglou, U.o.M. 2009
% contact at: vogia@yahoo.com
% ------------------------------------------------------------------------
    
if nargin == 3
    solver = 'fmincon';
end

if strcmp(solver,'fminunc')==1
    % the input vector theta should be unconstrained
    % but the log likelihood accepts only constrained parameters
    if strcmp(spec.purpose,'fitCopulaGARCH')==1
    spec.purpose = 'fitCopula';   
    theta = RescaleParameters(theta, 1, spec);
    spec.purpose = 'fitCopulaGARCH';
    else
    theta = RescaleParameters(theta, 1, spec);
    end
    % thus the initial theta is transformed to the corresponding
    % constrained theta
end
%theta % you might want to see this if the procedure crashes...
if max(max(data))>1 || min(min(data))<0
    error('data fed to copula should be uniform')
end
corrspec = spec.depspec;
% -----------------------------------------------------------------------
% the Clayton Copula Log Likelihood
% -----------------------------------------------------------------------
if strcmp(spec.family,'Clayton')==1
    T=size(data,1);
    u=data(:,1);
    v=data(:,2);
    if strcmp(corrspec,'static')==1 
    tau=theta; d=2*tau/(1-tau); d=repmat(d,[T,1]); 
    elseif strcmp(corrspec,'Patton')==1 
    tau=Pattoneq(theta,data,'Clayton'); d=2*tau./(1-tau);
    end
    out1=(u.*v).^(-1-d);
    out2=(u.^(-d)+v.^(-d)-1).^((-1./d) -2);
    out=(1+d).*out1.*out2;
    LL = -log(out);
    LogL=sum(-log(out));
    Rt=tau;
    varargout{1} = Rt;
end
% -----------------------------------------------------------------------
% the SJC copula log likelihood
% -----------------------------------------------------------------------
if strcmp(spec.family,'SJC')==1
T=size(data,1);
u=data(:,1); v=data(:,2);
if strcmp(corrspec,'static')==1 
    if size(theta,1)==2
    tauU=repmat(theta(1),[T,1]); tauL=repmat(theta(2),[T,1]);
    elseif size(theta,2)==2
    tauU=repmat(theta(1,1),[T,1]); tauL=repmat(theta(1,2),[T,1]);
    end
elseif strcmp(corrspec,'Patton')==1 
    tau=Pattoneq(theta,data,'SJC');
    tauU=tau(:,1);
    tauL=tau(:,2); 
end   
k1 =  1./log2(2-tauU);
k2 = -1./log2(tauL);
JC1=(k1.*k2.*(1 - 1./(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2)).^(1./k1 - 1).*(1./k2 + 1).*(1 - u).^(k1 - 1).*(1 - v).^(k1 - 1))./((1 - (1 - u).^k1).^(k2 + 1).*(1 - (1 - v).^k1).^(k2 + 1).*(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2 + 2));
JC2=(k1.*(1 - 1./(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2)).^(1./k1 - 2).*(1./k1 - 1).*(1 - u).^(k1 - 1).*(1 - v).^(k1 - 1))./((1 - (1 - u).^k1).^(k2 + 1).*(1 - (1 - v).^k1).^(k2 + 1).*(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(2./k2 + 2));
out1=JC1-JC2;

k1 =  1./log2(2-tauL);
k2 = -1./log2(tauU);
u  = 1-u;
v  = 1-v;
JC3=(k1.*k2.*(1 - 1./(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2)).^(1./k1 - 1).*(1./k2 + 1).*(1 - u).^(k1 - 1).*(1 - v).^(k1 - 1))./((1 - (1 - u).^k1).^(k2 + 1).*(1 - (1 - v).^k1).^(k2 + 1).*(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2 + 2));
JC4=(k1.*(1 - 1./(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(1./k2)).^(1./k1 - 2).*(1./k1 - 1).*(1 - u).^(k1 - 1).*(1 - v).^(k1 - 1))./((1 - (1 - u).^k1).^(k2 + 1).*(1 - (1 - v).^k1).^(k2 + 1).*(1./(1 - (1 - u).^k1).^k2 + 1./(1 - (1 - v).^k1).^k2 - 1).^(2./k2 + 2));
out2=JC3-JC4;
LL = - log(.5*(out1+out2));
LogL=-sum(log(.5*(out1+out2)));
Rt=[tauU tauL];
varargout{1} = Rt;
end
[T,N]=size(data);
% -----------------------------------------------------------------------
% the t - Copula log - likelihood
% -----------------------------------------------------------------------
if strcmp(spec.family,'t')==1
    nu=theta(1);
   
    trdata=tinv(data,nu);
    
    if strcmp(corrspec,'DCC')==1 
        [Rt, veclRt]=DCCeq(theta(2:end),trdata);
    end
    
    if strcmp(corrspec,'static')==1
    rrho=corr(data,'type','kendall'); %since data is U(0,1) corr(data) is the rank correlation
    rho=sin(.5*pi*rrho);
    Rt=repmat(rho,[1 1 T]);
    else
    varargout{1} = veclRt;    
    end
    varargout{2} = Rt;
    
    % The T Copula likelihood function 
    LL=zeros(T,1); 
    for i=1:T
        LL(i) = gammaln((nu+N)/2) + (N-1)*gammaln(nu/2) - N*gammaln((nu+1)/2) - 0.5*log(det(Rt(:,:,i)));
        LL(i) = LL(i) - (nu+N)/2*log(1+trdata(i,:)*inv(Rt(:,:,i))*trdata(i,:)'./nu);
        LL(i) = LL(i) + (nu+1)/2*sum(log(1+(trdata(i,:).^2/nu)));
    end
    LL = - LL;
    LogL=sum(LL);
end
% ------------------------------------------------------------------------
% the Gaussian Copula log likelihood
% ------------------------------------------------------------------------
if strcmp(spec.family,'Gaussian')==1 
    trdata=norminv(data);
    
    if strcmp(corrspec,'DCC')==1 
    [Rt, veclRt]=DCCeq(theta,trdata);
    elseif strcmp(corrspec,'TVC')==1
    [Rt, veclRt]=TVCeq(theta,trdata);
    end  
    varargout{1} = veclRt; varargout{2} = Rt;
    LL=zeros(T,1); 
    for i=1:T
        LL(i)=-.5*log(det(Rt(:,:,i)));
        LL(i)=LL(i)-.5*trdata(i,:)*(inv(Rt(:,:,i))-eye(N))*trdata(i,:)';    
    end
    LL = - LL;
    LogL = sum(LL);
end


  