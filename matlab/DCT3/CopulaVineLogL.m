function [LogL,LL]=CopulaVineLogL(phi,data,spec, solver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % the code segment for copula vine log - likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dec=spec.decomposition;  
[T,N]=size(data);
lik=zeros(N-1,N-1); ll = cell(N-1,N-1);
u=cell(N,N);
if min(min(data))<0 || max(max(data))>1
    data=empiricalCDF(data);
    display('The data you provided was not usiform, it was transformed to uniform with')
    display('the empiricalCDF.m function. Press alt ctrl to abord')
end
if strcmp(spec.family,'SJC')==0
    if size(phi,1)~=.5*N*(N-1) || size(phi,2)~=1
    error('the parameters vector should be a vector with .5*N*(N-1) rows')
    end
    %if strcmp(solver, 'fminunc') == 1
    % the theta vector that is fed to the log likelihood should contain
    % constrained parameters
    %phi = RescaleParams(phi, 1, spec);
    %end
    phi=crIPmat(phi);
if strcmp(dec,'CVine')==1
    for j=1:N-1
    for i=1:N-j
       if j==1
       [lik(j,i),ll{j,i}]=CopulaLogL(phi(j,i),[data(:,1),data(:,i+1)],spec,solver);
       else
       [lik(j,i),ll{j,i}]=CopulaLogL(phi(j,i),[u{j-1,1},u{j-1,i+1}],spec,solver);
       end
   end
   if j<N-1
       for i=1:N-j
       if j==1
       u{1,i}=hfunction(data(:,i+1),data(:,1),phi(1,i),spec,solver);
       else
       u{j,i}=hfunction(u{j-1,i+1},u{j-1,1},phi(j,i),spec,solver);
       end
       end
   end
   end
elseif strcmp(dec,'DVine')==1
    for i=1:N-1
    [lik(1,i), ll{1,i}]=CopulaLogL(phi(1,i),[data(:,i) data(:,i+1)],spec, solver);
    end
    u{1,1}=hfunction(data(:,1),data(:,2),phi(1,1),spec,solver);
    for k=1:N-3
    u{1,2*k}=hfunction(data(:,k+2),data(:,k+1),phi(1,k+1),spec,solver);
    u{1,2*k+1}=hfunction(data(:,k+1),data(:,k+2),phi(1,k+1),spec,solver);
    end
    u{1,2*N-4}=hfunction(data(:,N),data(:,N-1),phi(1,N-1),spec,solver);
    for j=2:(N-1)
    for i=1:(N-j)
        [lik(j,i),ll{j,i}]=CopulaLogL(phi(j,i),[u{j-1,2*i-1} u{j-1,2*i}],spec,solver);
    end
    if j<N-1
    u{j,1}=hfunction(u{j-1,1},u{j-1,2},phi(j,1),spec,solver); 
    if N>4
        for i=1:(N-j-2)
        u{j,2*i}=hfunction(u{j-1,2*i+2},u{j-1,2*i+1},phi(j,i+1),spec,solver);
        u{j,2*i+1}=hfunction(u{j-1,2*i+1},u{j-1,2*i+2},phi(j,i+1),spec,solver);
        end
    end
    u{j,2*N-2*j-2}=hfunction(u{j-1,2*N-2*j},u{j-1,2*N-2*j-1},phi(j,N-j),spec,solver);
    end
    end
end
LogL=sum(sum(lik));
else
% ----------------------  some comments ---------------------------------
% In the SJC copula decomposition, each bivariate copula has two
% parameters, unlike the t or Clayton copula which have only one. Therefore
% the output and initial values are not in a matrix, but in a cell array, with
% same dimension as the matrix. For example in the t - copula theta(1,1)
% contains the dof parameter for the first copula in the cascade, whereas
% in the SJC copula, theta{1,1} contains two values, the upper and lower
% tail dependence of the first copula.
% ------------------------------------------------------------------------
phi = reshape(phi,[.5*N*(N-1),2]);
phi=crIPmat(phi);
%xyz = icrIPmat(phi)
if strcmp(dec,'CVine')==1
    for j=1:N-1
    for i=1:N-j
       if j==1
       xy = phi{j,i};
       %if strcmp(solver,'fminunc')==1
       %xy = RescaleParams(xy,2,[],CopulaSpec);
       %end
       [lik(j,i),ll{j,i}]=CopulaLogL(xy,[data(:,1),data(:,i+1)],spec,solver);
       else
       xy = phi{j,i};
       %if strcmp(CopulaSpec.solver,'fminunc')==1
       %xy = RescaleParams(xy,2,[],CopulaSpec);
       %end
       [lik(j,i),ll{j,i}]=CopulaLogL(xy,[u{j-1,1},u{j-1,i+1}],spec,solver);
       end
   end
   if j<N-1
       for i=1:N-j
       if j==1
       yx = phi{1,i};
       %if strcmp(CopulaSpec.solver,'fminunc')==1
       %yx = RescaleParams(yx,2,[],CopulaSpec);
       %end
       u{1,i}=hfunction(data(:,i+1),data(:,1),yx,spec,solver);
       else
       yx = phi{j,i};
       %if strcmp(CopulaSpec.solver,'fminunc')==1
       %yx = RescaleParams(yx,2,[],CopulaSpec);
       %end
       u{j,i}=hfunction(u{j-1,i+1},u{j-1,1},yx,spec,solver);
       end
       end
   end
   end
elseif strcmp(dec,'DVine')==1
    for i=1:N-1
    xy = phi{1,i};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    [lik(1,i), ll{1,i}]=CopulaLogL(xy,[data(:,i) data(:,i+1)],spec,solver);
    end
    xy = phi{1,1};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    u{1,1}=hfunction(data(:,1),data(:,2),xy,spec,solver);
    for k=1:N-3
    xy = phi{1,k+1};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    u{1,2*k}=hfunction(data(:,k+2),data(:,k+1),xy,spec,solver);
    u{1,2*k+1}=hfunction(data(:,k+1),data(:,k+2),xy,spec,solver);
    end
    xy = phi{1,N-1};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    u{1,2*N-4}=hfunction(data(:,N),data(:,N-1),xy,spec,solver);
    for j=2:(N-1)
    for i=1:(N-j)
        xy = phi{j,i};
        %if strcmp(CopulaSpec.solver,'fminunc')==1
        %xy = RescaleParams(xy,2,[],CopulaSpec);
        %end
        [lik(j,i),ll{j,i}]=CopulaLogL(xy,[u{j-1,2*i-1} u{j-1,2*i}],spec,solver);
    end
    if j<N-1
    xy = phi{j,1};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    u{j,1}=hfunction(u{j-1,1},u{j-1,2},xy, spec, solver); 
    if N>4
        for i=1:(N-j-2)
        xy = phi{j,i+1};
        %if strcmp(CopulaSpec.solver,'fminunc')==1
        %xy = RescaleParams(xy,2,[],CopulaSpec);
        %end
        u{j,2*i}=hfunction(u{j-1,2*i+2},u{j-1,2*i+1},xy,spec, solver);
        u{j,2*i+1}=hfunction(u{j-1,2*i+1},u{j-1,2*i+2},xy,spec, solver);
        end
    end
    xy = phi{j,N-j};
    %if strcmp(CopulaSpec.solver,'fminunc')==1
    %    xy = RescaleParams(xy,2,[],CopulaSpec);
    %end
    u{j,2*N-2*j-2}=hfunction(u{j-1,2*N-2*j},u{j-1,2*N-2*j-1},xy,spec,solver);
    end
    end
end
LogL=sum(sum(lik));

end
LL = sum(icrIPmat(ll))';
end
  