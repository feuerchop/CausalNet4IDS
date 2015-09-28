function [A, B, lower, upper] = CreateFminconConstraints(spec)
% This function creates the arrays A, B, lower and upper that are used as
% inputs to fmincon.
purp = spec.purpose;

switch purp
    case 'fitGARCH'

        mp= spec.mparams; vp = spec.vparams; dp = spec.dparams;
        if dp ==  1 % upper and lower values for the dof parameter in tGARCH
            ldp = 2;
            udp = 200;
        elseif dp == 2.% upper and lower values for the dof and asymmetry parameter in SkewtGARCH
            ldp = [2.01;-.5];
            udp = [200; .5];
        else
            ldp = []; udp = [];
        end

        if vp == 3 % upper and lower values for the GARCH(1,1) model
            A = [zeros(1,mp) [0 1 1] zeros(1,dp)];
            lower = [-.5*ones(mp,1); [10^-3*spec.vtheta0(1) 10^-9 10^-9]';ldp];
            upper = [.5*ones(mp,1); [5*spec.vtheta0(1) .5 .9995]';udp];
            
        elseif vp == 4 % upper and lower values for the GJR(1,1) model
            A = [zeros(1,mp) [0 .5 1 .5] zeros(1,dp)];
            lower = [-.5*ones(mp,1); [10^-3*spec.vtheta0(1) 10^-9 10^-9 10^-9]';ldp];
            upper = [.5*ones(mp,1); [5*spec.vtheta0(1) .5 .9995 .5]';udp];
        end
        B = .99999;
    case 'fitCopula'
        if strcmp(spec.family,'t') == 1 && strcmp(spec.depspec,'static')==1
            A = 0; % the dof parameter of the static t copula, upper and lower values
            lower = 2.01;
            upper = 200;
        elseif strcmp(spec.family,'t') == 1 && strcmp(spec.depspec,'static')==0
            A = [0 1 1]; % upper and lower values for the tDCC copula parameters
            lower = [2.01; 10^-9; 10^-9];
            upper = [200; .5; .99];
        end
        if strcmp(spec.family,'Gaussian')==1
            A = [1 1]; % upper and lower values for the Gaussian DCC copula parameters
            lower = [10^-9; 10^-9];
            upper = [ .5; .9997];
        end
        if strcmp(spec.family,'Clayton') == 1 && strcmp(spec.depspec,'static')==1
            A = 0; % upper and lower value for the static Clayton parameter
            lower = 10^-12;
            upper = .85;
        elseif strcmp(spec.family,'Clayton') == 1 && strcmp(spec.depspec,'static')==0
            A = zeros(1,3); % upper and lower values for the tv Clayton copula parameters
            lower = -15*ones(3,1);
            upper = 15*ones(3,1);
        end
        if strcmp(spec.family,'SJC') == 1 && strcmp(spec.depspec,'static')==1
            A = zeros(1,2); % upper and lower values for the static SJC copula parameters
            lower = [10^-12; 10^-12];
            upper = [.85; .85];
        elseif strcmp(spec.family,'SJC') == 1 && strcmp(spec.depspec,'static')==0
            A = zeros(1,6); % upper and lower values for the tv SJC copula parameters
            lower = -10*ones(6,1);
            upper = 10*ones(6,1);
        end
        B = .9999;
    case 'fitCopulaGARCH'
        n = spec.size;
        m = spec.vecsize;
        mp= spec.mparams; vp = spec.vparams; dp = spec.dparams;
        if dp ==  1
            ldp = 2;
            udp = 200;
        elseif dp == 2
            ldp = [2.01;-.5];
            udp = [200; .5];
        else
            ldp = []; udp = [];
        end
        Am = zeros(n,n*m);
        if vp == 3
          
            for i=1:n
                Am(i,:) = [zeros(1,m*(i-1)) zeros(1,mp) [0 1 1] zeros(1,dp) zeros(1,n*m - i*m)];
            end
            
            lowerm = repmat([-.5*ones(mp,1); [10^-20 10^-9 10^-9]';ldp],[n,1]);
            upperm = repmat([.5*ones(mp,1); [10 .5 .9995]';udp],[n,1]);
            
        elseif vp == 4
            for i=1:n
                Am(i,:) = [zeros(1,m*(i-1)) zeros(1,mp) [0 .5 1 .5] zeros(1,dp) zeros(1,n*m - i*m)];
            end
            lowerm = repmat([-.5*ones(mp,1); [10^-20 10^-9 10^-9 10^-9]';ldp],[n,1]);
            upperm = repmat([.5*ones(mp,1); [10 .5 .9995 .5]';udp],[n,1]);
        end
        
        if strcmp(spec.family,'t') == 1 && strcmp(spec.depspec,'static')==1
            Ac = 0;
            lowerc = 2.01;
            upperc = 200;
        elseif strcmp(spec.family,'t') == 1 && strcmp(spec.depspec,'static')==0
            Ac = [0 1 1];
            lowerc = [2.01; 10^-9; 10^-9];
            upperc = [200; .08; .995];
        end
        if strcmp(spec.family,'Gaussian')==1
            Ac = [1 1];
            lowerc = [10^-9; 10^-9];
            upperc = [ .5; .9872];
        end
        if strcmp(spec.family,'Clayton') == 1 && strcmp(spec.depspec,'static')==1
            Ac = 0;
            lowerc = 10^-12;
            upperc = .85;
        elseif strcmp(spec.family,'Clayton') == 1 && strcmp(spec.depspec,'static')==0
            Ac = zeros(1,3); 
            lowerc = -10*ones(3,1);
            upperc = 10*ones(3,1);
        end
        if strcmp(spec.family,'SJC') == 1 && strcmp(spec.depspec,'static')==1
            Ac = zeros(1,2);
            lowerc = [10^-12; 10^-12];
            upperc = [.85; .85];
        elseif strcmp(spec.family,'SJC') == 1 && strcmp(spec.depspec,'static')==0
            Ac = zeros(1,6); 
            lowerc = -10*ones(6,1);
            upperc = 10*ones(6,1);
        end
        A = [Am zeros(n,size(Ac,2)); zeros(1,size(Am,2)) Ac];
        lower = [lowerm;lowerc];
        upper = [upperm;upperc];
        B = .99999*ones(n+1,1);
    case 'fitCopVine'
        n = spec.size;
        if strcmp(spec.family,'t')==1
            A = zeros(1,.5*n*(n-1));
            B = zeros;
            lower = 2.1*ones(.5*n*(n-1),1);
            upper = 200*ones(.5*n*(n-1),1);
        elseif strcmp(spec.family,'Clayton')==1
            A = [];
            B = [];
            lower = 10^-6*ones(.5*n*(n-1),1);
            upper = .85*ones(.5*n*(n-1),1);
        elseif strcmp(spec.family,'SJC')==1
            A = [];
            B = [];
            lower = 10^-6*ones(.5*n*(n-1),2);
            upper = .85*ones(.5*n*(n-1),2);
        end
            
end

        
        


           