function [out, firstDer, secondDer] = RescaleParameters(theta, flag, spec)

% transforms unconstrained parameters to constraint and vice versa
% USAGE: to estimate a copula GARCH model with fminunc;
% INPUTS:
% theta:        column vector of parameters from a supported model
% flag:         integer with values 1 or 2. Set 1 when you want to get the
%               actuall parameters from the unconstrained ones and 2 for
%               vice versa
% garchSpec:    structure that contains the model specifications 

% DETAILS:
% Let u be a vector of parameters obtained from fminunc. To get the actuall
% vector v, of parameters use:
% v = Rescaleparams(u,1,garchSpec,copulaspec)
% If you want to get u from v, use:
% u = Rescaleparams(v,2,garchSpec,copulaspec)
% to constrain the parameters the following functions are used:
% For params in (0,a) ---> v = a/(1+exp(-u))
% For params in (a, +inf) ----> v = a + exp(u);
% For params in (-a, a) ---> v = -a +2*a*exp(u)/(1+exp(u))
n = spec.size;
purp = spec.purpose;
%%%%%%%%%%%%%%%%%%%%%%% ---- Only GARCH ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(purp,'fitGARCH')==1 
    
    %ss = spec.vecsize;
%    theta = reshape(theta,[ss,n]); % input theta is a column vec
    mp = spec.mparams; vp = spec.vparams; 
    
    switch flag
        case 1 % create constrained parameters from unconstrained
            
            % separate theta to mean, variance and distributional params
            
            uncmparams = theta(1:mp,:);  
            uncvparams = theta(mp+1:mp+vp,:);   consvparams = zeros(size(uncvparams));
            uncdparams = theta(mp+vp+1:end);   consdparams = zeros(size(uncdparams));
            
            % mean parameters are actually unconstrained but I will keep
            % them in (-.5, .5) to avoid explosive behaviour from the
            % fminunc
            
            consmparams = FromUnc2Cons(uncmparams,.5,3);
            
            % varparams are [omega alfa beta] or [omega alfa beta gamma]
            % depending on spec.VarEq. The following hold
            % omega > 0
            % 1> alfa > 0
            % 1> beta > 0
            % gamma > -alfa for GJR
            % alfa + beta < 1
            % beta + .5*(alfa + gamma) < 1 for GJR
            
            consvparams(1) = FromUnc2Cons(uncvparams(1),0,2);
            consvparams(2) = FromUnc2Cons(uncvparams(2),.5,1); % alfa is kept in (0,.5)
            consvparams(3) = FromUnc2Cons(uncvparams(3),.9995,1); % beta is kept in (0, .9995)
            
            if strcmp(spec.VarEq,'GARCH(1,1)') == 1

                    if consvparams(2) + consvparams(3) > .99999;
                        consvparams(2) = .99999 - consvparams(3);
                    end

            elseif strcmp(spec.VarEq,'GJR(1,1)') == 1
                consvparams(4) = FromUnc2Cons(uncvparams(4),.25,3);

                    % check if gamma > - alfa, if not, set alfa = - gamma
                    if consvparams(4) < - consvparams(2)
                       consvparams(4) = - consvparams(2);
                    end
                    % check if beta + 2*(alfa + gamma) < 1
                    if consvparams(3)+ .5*(consvparams(2)+consvparams(4))>.99999
                        consvparams(2) = 2*(.99999 - consvparams(3)) - consvparams(4);
                    end

            end
            
            % for the distributional parameters, the following restrictions hold:
            % dof > 2
            % lambda in (-1,1)
            
            if strcmp(spec.distr,'Gaussian')==1
                consdparams = consvparams; % both are an empty matrix
            else
                consdparams(1) = FromUnc2Cons(uncdparams(1),[2 98],4);
                % the dof parameter is held in (2,100)
            end
            
            if strcmp(spec.distr,'SkewT')==1
                % again to avoid explosive behaviour from fminunc lambda is
                % kept in (-0.5, 0.5)
                consdparams(2,:) = FromUnc2Cons(uncdparams(2),.5,3);
            end
            
            out = [consmparams;consvparams;consdparams];
            firstDer = [];
            secondDer = [];
            
        case 2 %create unconstrained from constrained parameters
             
            % separate theta to mean, variance and distributional params
            
            consmparams = theta(1:mp,:);          
            consvparams = theta(mp+1:mp+vp,:);   uncvparams = zeros(size(consvparams)); d1vp = uncvparams; d2vp = d1vp;
            consdparams = theta(mp+vp+1:end);    uncdparams = zeros(size(consdparams)); d1dp = uncdparams; d2dp = d1dp;
            
            % the functions that create the unconstrained parameters from
            % the constrained ones are the inverse functions of the ones
            % used to get the constrained parameters, thus
            % x = log(y/(a-y)) belongs to R if  y belongs to (0,a)
            % x = log(y-a) belongs to R if y>a
            % x = log((a+y)/(a-y)) belongs to R if  -a < y < a
            % x = log((a-y)/(y-a-b)) belongs in R if a < y < a+b
            
            [uncmparams, d1mp, d2mp] = FromCons2Unc(consmparams,.5,3);
            
            % unconstrained variance parameters
            
            [uncvparams(1,:), d1vp(1,:), d2vp(1,:)] = FromCons2Unc(consvparams(1,:),0,2);
            [uncvparams(2,:), d1vp(2,:), d2vp(2,:)] = FromCons2Unc(consvparams(2,:),0.5,1);
            [uncvparams(3,:), d1vp(3,:), d2vp(3,:)] = FromCons2Unc(consvparams(3,:),0.9995,1);
            
            if strcmp(spec.VarEq,'GJR(1,1)')==1
                
            [uncvparams(4,:), d1vp(4,:), d2vp(4,:)] = FromCons2Unc(consvparams(4,:),0.25,3);
            
            end
            % unconstrained  distr. parameters
            
            if strcmp(spec.distr,'Gaussian')==1
                uncdparams = consdparams;
            else
                [uncdparams(1,:), d1dp(1,:), d2dp(1,:)] = FromCons2Unc(consdparams(1,:),[2 98],4);
            end
            
            if strcmp(spec.distr,'SkewT')==1
                [uncdparams(2,:), d1dp(2,:), d2dp(2,:)] = FromCons2Unc(consdparams(2,:),.5,3);
            end
            % concatenate results to create the output
            out = [uncmparams;uncvparams;uncdparams];
            firstDer = [d1mp; d1vp; d1dp];
            secondDer = [d2mp; d2vp; d2dp];
    end
%%%%%%%%%%%%% --------- Only Copula ----------- %%%%%%%%%%%%%%%%%%%%%
elseif strcmp(purp,'fitCopula')==1 || strcmp(purp,'fitCopVine')==1
    
    switch flag
        case 1 % create constrained parameters
            unccparams = theta; conscparams = zeros(size(unccparams));
            % t copula
            if strcmp(spec.family,'t')==1 % the t copula parameters
                if strcmp(spec.depspec,'static')==1
                conscparams = FromUnc2Cons(unccparams,[2 198],4);
                elseif strcmp(spec.depspec , 'DCC')==1
                    conscparams(1,:) = FromUnc2Cons(unccparams(1,:),[2 198],4);
                    conscparams(2,:) = FromUnc2Cons(unccparams(2,:),.5,1);
                    conscparams(3,:) = FromUnc2Cons(unccparams(3,:),0.999,1);
                    if conscparams(2) + conscparams(3) > .99999;
                        conscparams(2) = .99999 - conscparams(3);
                    end
                end
            end
            if strcmp(spec.family,'Gaussian')==1
            % gaussian copula
            conscparams(1,:) = FromUnc2Cons(unccparams(1,:),0.5,1);
            conscparams(2,:) = FromUnc2Cons(unccparams(2,:),0.999,1);
            if conscparams(1) + conscparams(2) > .99999;
                conscparams(1) = .99999 - conscparams(2);
            end
            end
            
            % tv SJC and Claton have unconstrained parameters 
            
            if (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==1
                
                conscparams = FromUnc2Cons(unccparams,.85,1); % keep the param in (0, 0.85)
            elseif (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==0
                conscparams = FromUnc2Cons(unccparams,10,3);
            end
            out = conscparams; firstDer = []; secondDer = [];
            
        case 2 % create the unconstrained parameters
            
            conscparams = theta; unccparams = zeros(size(conscparams)); d1cp = unccparams; d2cp = d1cp;
            
            % t copula
            if strcmp(spec.family,'t')==1 % the t copula parameters
                if strcmp(spec.depspec,'static')==1
                    [unccparams, d1cp, d2cp] = FromCons2Unc(conscparams,[2,198],4);
                elseif strcmp(spec.depspec , 'DCC')==1
                    [unccparams(1,:), d1cp(1,:), d2cp(1,:)] = FromCons2Unc(conscparams(1,:),[2,198],4);
                    [unccparams(2,:), d1cp(2,:), d2cp(2,:)] = FromCons2Unc(conscparams(2,:),0.5,1);
                    [unccparams(3,:), d1cp(3,:), d2cp(3,:)] = FromCons2Unc(conscparams(3,:),0.999,1);
                end
            end
            if strcmp(spec.family,'Gaussian')==1
            % gaussian copula
            [unccparams(1,:), d1cp(1,:), d2cp(1,:)] = FromCons2Unc(conscparams(1,:),0.5,1);
            [unccparams(2,:), d1cp(2,:), d2cp(2,:)] = FromCons2Unc(conscparams(2,:),0.999,1);
            end
            
            % tv SJC and Claton have unconstrained parameters 
            
            if (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==1
                
                [unccparams, d1cp, d2cp] = FromCons2Unc(conscparams,.85,1); % keep the param in (0, 0.85)
            elseif  (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==0
                [unccparams, d1cp, d2cp] = FromCons2Unc(conscparams,10,3);
            end
            
            out = unccparams;
            firstDer = d1cp; secondDer = d2cp;
    end
%%%%%%%%%%%%%%%%%%%% ----- Copula GARCH one step --------- %%%%%%%%%%%%%%%%

elseif strcmp(purp,'fitCopulaGARCH')==1
    
    ss = spec.vecsize;
    phi = theta(n*ss+1:end); % copula params
    theta = theta(1:n*ss); % GARCHparams  
    theta = reshape(theta,[ss,n]); % input theta is a column vec
    mp = spec.mparams; vp = spec.vparams;
    
     switch flag
        case 1 % create constrained parameters from unconstrained
            
            % separate theta to mean, variance and distributional params
            
            uncmparams = theta(1:mp,:);  
            uncvparams = theta(mp+1:mp+vp,:);   consvparams = zeros(size(uncvparams));
            uncdparams = theta(mp+vp+1:end,:);   consdparams = zeros(size(uncdparams));
            
            % mean parameters are actually unconstrained but I will keep
            % them in (-.5, .5) to avoid explosive behaviour from the
            % fminunc
            
            consmparams = FromUnc2Cons(uncmparams,.5,3);
            
            % varparams are [omega alfa beta] or [omega alfa beta gamma]
            % depending on spec.VarEq. The following hold
            % omega > 0
            % 1> alfa > 0
            % 1> beta > 0
            % gamma > -alfa for GJR
            % alfa + beta < 1
            % 2*(alfa + gamma) + beta < 1 for GJR
            
            consvparams(1,:) = FromUnc2Cons(uncvparams(1,:),0,2);
            consvparams(2,:) = FromUnc2Cons(uncvparams(2,:),.5,1); % alfa is kept in (0,.5)
            consvparams(3,:) = FromUnc2Cons(uncvparams(3,:),.9995,1); % beta is kept in (0, .9995)
            
            if strcmp(spec.VarEq,'GARCH(1,1)') == 1
                for j = 1:n
                    if consvparams(2,j) + consvparams(3,j) > .99999;
                        consvparams(2,j) = .99999 - consvparams(3,j);
                    end
                end
            elseif strcmp(spec.VarEq,'GJR(1,1)') == 1
                consvparams(4,:) =FromUnc2Cons(uncvparams(4,:),.25,3); % kept in (-.25,.25)
                for j=1:n
                    % check if gamma > - alfa, if not, set alfa = - gamma
                    if consvparams(4,j) < - consvparams(2,j)
                       consvparams(4,j) =  -consvparams(2,j);
                    end
                    % check if beta + .5*(alfa + gamma) < 1
                    if consvparams(3,j)+ .5*(consvparams(2,j)+consvparams(4,j))>.99999
                        consvparams(2,j) = 2*(.99999 - consvparams(3,j)) - consvparams(4,j);
                    end
                end
            end
            
            % for the distributional parameters, the following restrictions hold:
            % dof > 2
            % lambda in (-1,1)
            
            if strcmp(spec.distr,'Gaussian')==1
                consdparams = consvparams; % both are an empty matrix
            else
                consdparams(1,:) = FromUnc2Cons(uncdparams(1,:),[2 100],4);
            end
            
            if strcmp(spec.distr,'SkewT')==1
                % again to avoid explosive behaviour from fminunc lambda is
                % kept in (-0.5, 0.5)
                consdparams(2,:) = FromUnc2Cons(uncdparams(2,:),.5,3);
            end
            
            %%%%%%%%%%% --------  Copula parameters -----------------------
            
            unccparams = phi; conscparams = zeros(size(unccparams));
            % t copula
            if strcmp(spec.family,'t')==1 % the t copula parameters
                conscparams(1,:) = FromUnc2Cons(unccparams(1,:),[2,198],4);
                if strcmp(spec.depspec , 'DCC')==1
                    conscparams(2,:) = FromUnc2Cons(unccparams(2,:),0.5,1);
                    conscparams(3,:) = FromUnc2Cons(unccparams(3,:),0.999,1);
                    if conscparams(2) + conscparams(3)>.99999;
                        conscparams(2) = .99999 - conscparams(3);
                    end
                end
            end
            % gaussian copula
            if strcmp(spec.family,'Gaussian')==1 % the t copula parameters
            conscparams(1,:) = FromUnc2Cons(unccparams(1,:),0.5,1);
            conscparams(2,:) = FromUnc2Cons(unccparams(2,:),0.999,1);
                 if conscparams(1) + conscparams(2)>.99999;
                        conscparams(1) = .99999 - conscparams(2);
                 end
            end
            
            % tv SJC and Claton have unconstrained parameters 
            
            if (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==1
                
                conscparams = FromUnc2Cons(unccparams,.85,1); % keep the param in (0, 0.85)
            elseif (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==0
                conscparams = FromUnc2Cons(unccparams,10,3);
                
            end
            
            out = reshape([consmparams;consvparams;consdparams],[n*ss,1]);
            out = [out;conscparams];
            firstDer = [];
            secondDer = [];
            
        case 2
            consmparams = theta(1:mp,:);          
            consvparams = theta(mp+1:mp+vp,:);   uncvparams = zeros(size(consvparams)); d1vp = uncvparams; d2vp = d1vp;
            consdparams = theta(mp+vp+1:end,:);    uncdparams = zeros(size(consdparams)); d1dp = uncdparams; d2dp = d1dp;
            conscparams = phi;     unccparams = zeros(size(conscparams)); d1cp = unccparams; d2cp = d1cp;
            
            % the functions that create the unconstrained parameters from
            % the constrained ones are the inverse functions of the ones
            % used to get the constrained parameters, thus
            % x = log(y/(a-y)) belongs to R if  y belongs to (0,a)
            % x = log(y-a) belongs to R if y>a
            % x = log((a+y)/(a-y)) belongs to R if  -a < y < a
            
            [uncmparams, d1mp, d2mp] = FromCons2Unc(consmparams,.5,3);
            
            % unconstrained variance parameters
            
            [uncvparams(1,:), d1vp(1,:), d2vp(1,:)] = FromCons2Unc(consvparams(1,:),0,2);
            [uncvparams(2,:), d1vp(2,:), d2vp(2,:)] = FromCons2Unc(consvparams(2,:),0.5,1);
            [uncvparams(3,:), d1vp(3,:), d2vp(3,:)] = FromCons2Unc(consvparams(3,:),0.9995,1);
            
            if strcmp(spec.VarEq,'GJR(1,1)')==1
                
            [uncvparams(4,:), d1vp(4,:), d2vp(4,:)] = FromCons2Unc(consvparams(4,:),0.25,3);
            
            end
            % unconstrained  distr. parameters
            
            if strcmp(spec.distr,'Gaussian')==1
                uncdparams = consdparams;
            else
                [uncdparams(1,:), d1dp(1,:), d2dp(1,:)] = FromCons2Unc(consdparams(1,:),[2 100],4);
            end
            
            if strcmp(spec.distr,'SkewT')==1
                [uncdparams(2,:), d1dp(2,:), d2dp(2,:)] = FromCons2Unc(consdparams(2,:),.5,3);
            end
            
            %%%% ------------ copula parameters --------------------
            
            conscparams = phi; 
            
            % t copula
            if strcmp(spec.family,'t')==1 % the t copula parameters
                [unccparams(1,:), d1cp(1,:), d2cp(1,:)] = FromCons2Unc(conscparams(1,:),[2 198],4);
                if strcmp(spec.depspec , 'DCC')==1
                    [unccparams(2,:), d1cp(2,:), d2cp(2,:)] = FromCons2Unc(conscparams(2,:),0.5,1);
                    [unccparams(3,:), d1cp(3,:), d2cp(3,:)] = FromCons2Unc(conscparams(3,:),0.999,1);
                end
            end
            % gaussian copula
            if strcmp(spec.family,'Gaussian')==1
            [unccparams(1,:), d1cp(1,:), d2cp(1,:)] = FromCons2Unc(conscparams(1,:),0.5,1);
            [unccparams(2,:), d1cp(2,:), d2cp(2,:)] = FromCons2Unc(conscparams(2,:),0.999,1);
            end
            % tv SJC and Claton have unconstrained parameters 
            
            if (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==1
                
                [unccparams, d1cp, d2cp] = FromCons2Unc(conscparams,.85,1); % keep the param in (0, 0.85)
                 elseif  (strcmp(spec.family,'Clayton')==1 || strcmp(spec.family,'SJC')==1) && strcmp(spec.depspec , 'static')==0
                [unccparams, d1cp, d2cp] = FromCons2Unc(conscparams,10,3);
                
            end
            
            out = reshape([uncmparams;uncvparams;uncdparams],[n*ss,1]);
            out = [out; unccparams];
            firstDer = [reshape([d1mp;d1vp;d1dp],[n*ss,1]);d1cp];
            secondDer = [reshape([d2mp;d2vp;d2dp],[n*ss,1]);d2cp];
     end
%elseif strcmp(purp,'fitCopVine')==1
    
end
            
                
            
            
            
            
            
            
            
