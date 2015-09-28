function tau=Pattoneq(theta,data,type)
T = size(data,1);
x = data(:,1);
y = data(:,2);
if  strcmp(type,'SJC')==1 
if size(theta,1) == 6
    theta = reshape(theta,[3,2]);
end
tau = ones(T,2); psi=ones(T,2);
tau(1,:) = .15*ones(1,2);
for jj = 2:T
    if jj<=10
        psi(:,1) = theta(1,1) + theta(2,1)*mean(abs(x(1:jj-1)-y(1:jj-1))) + theta(3,1)*psi(jj-1,1);
        psi(:,2) = theta(1,2) + theta(2,2)*mean(abs(x(1:jj-1)-y(1:jj-1))) + theta(3,2)*psi(jj-1,2);
    else
        psi(:,1) = theta(1,1) + theta(2,1)*mean(abs(x(jj-10:jj-1)-y(jj-10:jj-1))) + theta(3,1)*psi(jj-1,1);
        psi(:,2) = theta(1,2) + theta(2,2)*mean(abs(x(jj-10:jj-1)-y(jj-10:jj-1))) + theta(3,2)*psi(jj-1,2);
    end
   tau(jj,:) =.001+.75./(1+exp(-psi(jj,:)));	
end
elseif strcmp(type,'Clayton')==1
tau = ones(T,1); psi=zeros(T,1); tt = corr(data,'type','Kendall');
tau(1) = tt(1,2);
for jj = 2:T
   psi(jj) = theta(1,1) + theta(2,1)*abs(x(jj-1)-y(jj-1)) + theta(3,1)*psi(jj-1,1);
   tau(jj) =.0001+.75./(1+exp(-psi(jj)));
end

end