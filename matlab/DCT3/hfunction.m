function out=hfunction(u,v,theta,spec,solver)
% this function calculates the h function of some supported Copulas as
% introduced in Aas et al:"Pair - copula construction of multiple
% dependence"
% INPUTS:
% u,v:          uniform variables
% theta:        the copula parameters

[R,C]=size(theta);
if C>2
    error('theta is a scalar or a matrix with two columns at most')
end
if strcmp(solver,'fminunc')==1
theta = RescaleParameters(theta, 1, spec);
end

type=spec.family;
if strcmp(type, 't')==1
    if isscalar(theta)==1
    nu=theta;
    tau=corr(u,v,'type','kendall');
    rho=sin(tau*pi/2);
    elseif size(theta,1)==2
    nu=theta(1); rho=theta(2);
    elseif size(theta,1)==3
    nu=theta(1); 
    if strcmp(CopulaSpec.depspec,'DCC')==1
        [Rt, veclRt]=DCCeq(theta(2:end),[u,v],CopulaSpec.optimizer);
        rho=veclRt;
    elseif strcmp(CopulaSpec.depspec,'TVC')==1
        [Rt, veclRt]=TVCeq(theta(2:end),[u,v],CopulaSpec.optimizer);
        rho=veclRt;
    end
    end
%    if strcmp(CopulaSpec.solver,'fminunc')==1 && strcmp(CopulaSpec.use,'EstimateStartinVals')==0
%        nu=2.01+exp(nu); 
%    end
    out1=tinv(u,nu)-rho.*tinv(v,nu);
    out2=sqrt(((nu+tinv(v,nu).^2).*(1-rho.^2))./(nu+1));
    out=tcdf(out1./out2,nu+1);
elseif strcmp(type, 'Clayton')==1
    %if strcmp(CopulaSpec.solver,'fminunc')==1&& strcmp(CopulaSpec.use,'EstimateStartinVals')==0
    %    tau=.85./(1+exp(-theta));
    %else
        tau=theta;
    %end
    d=2*tau./(1-tau);
    out1=v.^(-d-1);
    out2=(u.^(-d)+v.^(-d)-1).^(-1-1./d);
    out=out1.*out2;
elseif strcmp(type, 'SJC')==1
    %if strcmp(CopulaSpec.solver,'fminunc')==1 && strcmp(CopulaSpec.use,'EstimateStartinVals')==0
    %    theta=.85./(1+exp(-theta)); 
    %end
    out1=hfuncJC(u,v,theta);
    out2=hfuncJC(1-u,1-v,[theta(2);theta(1)]);
    out=.5*(out1 - out2 + 1);
end
%clear rounding erros
T=size(out,1);
for i=1:T
    if out(i)>.9999
        out(i)=.9999;
    end
end
