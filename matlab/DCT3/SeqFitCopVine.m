function [theta, LogL]=SeqFitCopVine(data,spec)
%%%-----------Estimation of copula vine parameter in many steps-----------
% This function estimates the parameters of a Copula vine by using a
% multistep procedure. Each bivariate copula in the cascade is estimated
% sequentially
% INPUTS:
% data:             A TxN matrix of U(0,1) 
% spec:             Structure array with optimization parameters. To obtain
%                   first run the setCopulaVineLLinputs.m
% OUTPUT
% logL:             The log likelihood at the optimum
% theta:            The estimated parameters in vector form. 
% uncparams:        The unconstrained parameters
% WARNING:          This model unlike the others is estimated with fmincon,
%                   only. That is because the theory behind standard errors of
%                   sequentially estimated copula vines, is not established
%                   yet. Use the output of this as initial values for
%                   fitModel.m with spec.purpose = 'fitCopVine'
% ------------------------------------------------------------------------
% Author: Manthos Vogiatzoglou, UoM, 2008 - 2009
% contact at: vogia@yahoo.com
% -------------------------------------------------------------------------
display('**************************************************************')
display('the primary use of this function is to calculate the initial values for a Vine')
display('---------------------------------------------------------------------------------')
pause(2)
display('however the fitting a vine is very slow if t - copulas are assumed')
display('and the values obtained from the two methods are usually very very close')
display('----------------------------------------------------------------------')
pause(2)
display('therefore if you favor speed over accuracy, you may use this output')
display('----------------------------------------------------------------------')

% error checking and array initialization
if min(min(data))<0 || max(max(data))>1
    data=empiricalCDF(data);
    display('your data is transformed to uniform with the empiricalCDF function')
end
tic
N=size(data,2);
lik=zeros(N-1,N-1);
u=cell(N,N);
theta=zeros(N-1,N-1);
% -------------- the canonical vine decomposition ------------------------
if strcmp(spec.decomposition,'CVine')==1 && strcmp(spec.family,'SJC')==0

for j=1:N-1
   fprintf(1,'computing copula parameters for level %d\n',j)
   for i=1:N-j
      if j==1
      [theta(j,i),lik(j,i)]=fitModel(spec,[data(:,1),data(:,i+1)],'fmincon');
      else
      [theta(j,i),lik(j,i)]=fitModel(spec,[u{j-1,1}, u{j-1,i+1}],'fmincon');
      end
   end
   
   if j<N-1
       for i=1:N-j
       if j==1
       u{j,i}=hfunction(data(:,i+1),data(:,1),theta(j,i),spec,'fmincon');    
       else    
       u{j,i}=hfunction(u{j-1,i+1},u{j-1,1},theta(j,i),spec,'fmincon');
       end
       end
   end
end
theta=icrIPmat(theta);
% ---------------- canonical vine - SJC copula ----------------------------
elseif strcmp(spec.decomposition,'CanonicalVine')==1 && strcmp(spec.family,'SJC')==1
theta=cell(N-1,N-1);

for j=1:N-1
   fprintf(1,'computing copula parameters for level %d\n',j)
   for i=1:N-j
      if j==1
      [theta{j,i},lik(j,i)]=fitModel(spec,[data(:,1),data(:,i+1)],'fmincon');
      else
      [theta{j,i},lik(j,i)]=fitModel(spec,[u{j-1,1},u{j-1,i+1}],'fmincon');
      end
   end
   
   if j<N-1
       for i=1:N-j
       if j==1
       u{j,i}=hfunction(data(:,i+1),data(:,1),theta{j,i},spec,'fmincon');    
       else    
       u{j,i}=hfunction(u{j-1,i+1},u{j-1,1},theta{j,i},spec,'fmincon');
       end
       end
   end
end
theta=icrIPmat(theta);
% -----------------------------------------------------------------------
% the multi - step code for the D - Vine structure
% -----------------------------------------------------------------------
elseif strcmp(spec.decomposition,'DVine')==1 && strcmp(spec.family,'SJC')==0
    display('computing copula parameters for level 1');
for i=1:N-1
    [theta(1,i),lik(1,i)]=fitModel(spec, [data(:,i),data(:,i+1)], 'fmincon');
end
u{1,1}=hfunction(data(:,1),data(:,2),theta(1,1),spec,'fmincon');
for k=1:N-3
    u{1,2*k}=hfunction(data(:,k+2),data(:,k+1),theta(1,k+1),spec,'fmincon');
    u{1,2*k+1}=hfunction(data(:,k+1),data(:,k+2),theta(1,k+1),spec,'fmincon');
end
u{1,2*N-4}=hfunction(data(:,N),data(:,N-1),theta(1,N-1),spec,'fmincon');
for j=2:(N-1)
    fprintf(1,'computing copula parameters for level %d\n',j)
    for i=1:(N-j)
        [theta(j,i),lik(j,i)]=fitModel(spec, [u{j-1,2*i-1} u{j-1,2*i}],'fmincon');
    end
    if j<N-1
    u{j,1}=hfunction(u{j-1,1},u{j-1,2},theta(j,1),spec,'fmincon'); 
    if N>4
        for i=1:(N-j-2)
        u{j,2*i}=hfunction(u{j-1,2*i+2},u{j-1,2*i+1},theta(j,i+1),spec,'fmincon');
        u{j,2*i+1}=hfunction(u{j-1,2*i+1},u{j-1,2*i+2},theta(j,i+1),spec,'fmincon');
        end
    end
    u{j,2*N-2*j-2}=hfunction(u{j-1,2*N-2*j},u{j-1,2*N-2*j-1},theta(j,N-j),spec,'fmincon');
    end
end
theta=icrIPmat(theta);
% -------------------- D Vine SJC decomposition ---------------------------
elseif strcmp(spec.decomposition,'DVine')==1 && strcmp(spec.family,'SJC')==1
display('computing copula parameters for level 1');
theta=cell(N-1,N-1);

for i=1:N-1
    [theta{1,i},lik(1,i)]=fitModel(spec,[data(:,i),data(:,i+1)],'fmincon');
end
u{1,1}=hfunction(data(:,1),data(:,2),theta{1,1},spec,'fmincon');
for k=1:N-3
    u{1,2*k}=hfunction(data(:,k+2),data(:,k+1),theta{1,k+1},spec,'fmincon');
    u{1,2*k+1}=hfunction(data(:,k+1),data(:,k+2),theta{1,k+1},spec,'fmincon');
end
u{1,2*N-4}=hfunction(data(:,N),data(:,N-1),theta{1,N-1},spec,'fmincon');
for j=2:(N-1)
    fprintf(1,'computing copula parameters for level %d\n',j)
    for i=1:(N-j)
        [theta{j,i},lik(j,i)]=fitModel(spec,[u{j-1,2*i-1} u{j-1,2*i}],'fmincon');
    end
    if j<N-1
    u{j,1}=hfunction(u{j-1,1},u{j-1,2},theta{j,1},spec,'fmincon'); 
    if N>4
        for i=1:(N-j-2)
        u{j,2*i}=hfunction(u{j-1,2*i+2},u{j-1,2*i+1},theta{j,i+1},spec,'fmincon');
        u{j,2*i+1}=hfunction(u{j-1,2*i+1},u{j-1,2*i+2},theta{j,i+1},spec,'fmincon');
        end
    end
    u{j,2*N-2*j-2}=hfunction(u{j-1,2*N-2*j},u{j-1,2*N-2*j-1},theta{j,N-j},spec,'fmincon');
    end
end
theta=icrIPmat(theta);

end
LogL = sum(sum(lik));
%uncparams = [];
toc