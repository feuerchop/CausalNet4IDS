function [Rt, veclRt]=DCCeq(theta,data)
% Calculate the time varying correlations based on the DCC(1,1)
% specification of Engle, given a parameters vector theta and a matrix of
% standardized iid residuals data

[T,N]=size(data);
a=theta(1); b=theta(2);
if min(a,b)<0 || max(a,b)>1 || a+b > .999999;
    a = .9999 - b;
end
if max(max(isnan(data)))==1 || max(max(isinf(data)))==1
    error('data is NaN or inf')
end
%if strcmp(optimizer,'fmincon')==1
%    a=theta(1); b=theta(2);
%elseif strcmp(optimizer,'fminunc')==1
%    a=theta(1)^2/(1+theta(1)^2+theta(2)^2);
%    b=theta(2)^2/(1+theta(1)^2+theta(2)^2);
%end
Qt=zeros(N,N,T);
Qt(:,:,1)=cov(data);
if isnan(Qt(:,:,1))==1
    error('Q1 is nan')
end
Rt=zeros(N,N,T);
veclRt=zeros(T,N*(N-1)/2);
Rt(:,:,1)=corr(data);
for j=2:T
   Qt(:,:,j)=Qt(:,:,1)*(1-a-b);
   Qt(:,:,j)=Qt(:,:,j)+a*(data(j-1,:)'*data(j-1,:));
   Qt(:,:,j)=Qt(:,:,j)+b*Qt(:,:,j-1);
   Rt(:,:,j)=Qt(:,:,j)./(sqrt(diag(Qt(:,:,j)))*sqrt(diag(Qt(:,:,j)))');   
end
for j=1:T
veclRt(j,:)=vecl(Rt(:,:,j))';
if isnan(veclRt(j,:))==1
    Rt(:,:,j)=[1 .6;.6 1];
    cov(data)
end
end
