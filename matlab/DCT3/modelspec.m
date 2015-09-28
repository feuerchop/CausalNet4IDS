function spec = modelspec(data)
% creates the structure that defines the model characteristics
[T,n] = size(data);
spec.size = n;
if n == 1
    spec.purpose = 'fitGARCH';
        ab=input('input the lag - length of the AR terms in the mean equation and press enter:');
        spec.mparams=ab+1;
    	spec.mtheta0=[mean(data(:,1));zeros(ab,1)];
        ac=menu('define the variance equation','GARCH(1,1)','GJR(1,1)');
        if ac==1
            spec.VarEq='GARCH(1,1)';
            spec.vparams=3;
            spec.vtheta0=[.02*var(data(:,1));.08;.91];
        elseif ac==2
            spec.VarEq='GJR(1,1)';
            spec.vparams=4;
            spec.vtheta0=[.02*var(data(:,1));.05;.85;.15];
        end
        ad=menu('define the distribution of the residuals','Gaussian','T','skewT');
        if ad==1
            spec.distr='Gaussian';
            spec.dparams=0;
            spec.dtheta0=[];
        elseif ad==2
            spec.distr='T';
            spec.dparams=1;
            spec.dtheta0=10;
        else
            spec.distr='SkewT';
            spec.dparams=2;
            spec.dtheta0=[7;0];  
        end
        spec.vecsize = spec.mparams + spec.vparams + spec.dparams;
        ae=menu('define the PIT method','IFM','CML');
        if ae == 1
            spec.PIT = 'IFM';
        else
            spec.PIT = 'CML';
        end
end

if n > 1
    
    aa = menu('what do you want to estimate?','GARCH model for each series','Copula','Copula in one step','Copula Vine','Copula Vine sequentially');
    if aa == 1
        spec.purpose = 'fitGARCH';
    elseif aa == 2
        spec.purpose = 'fitCopula';
        spec.warning = 'this assumes that you have uniform data';
    elseif aa == 3
        spec.purpose = 'fitCopulaGARCH';
        spec.warning = 'this is for one step estimation procedure';
    elseif aa ==4
        if n > 2
            spec.purpose = 'fitCopVine';
        else
            error('Copula vines are not suitable for bivariate data')
        end
    elseif aa == 5
        spec.purpose = 'fitCopula';
        spec.comment = 'fit Copula Vine sequentially';
    end
    
    if aa==1 || aa==3  % when GARCH models will be estimated

        ab=input('input the lag - length of the AR terms in the mean equation and press enter:');
        spec.mparams=ab+1;
    	spec.mtheta0=[mean(data(:,1));zeros(ab,1)];
        ac=menu('define the variance equation','GARCH(1,1)','GJR(1,1)');
        if ac==1
            spec.VarEq='GARCH(1,1)';
            spec.vparams=3;
            spec.vtheta0=[.02*var(data(:,1));.08;.91];
        elseif ac==2
            spec.VarEq='GJR(1,1)';
            spec.vparams=4;
            spec.vtheta0=[.02*var(data(:,1));.05;.85;.15];
        end
        ad=menu('define the distribution of the residuals','Gaussian','T','skewT');
        if ad==1
            spec.distr='Gaussian';
            spec.dparams=0;
            spec.dtheta0=[];
        elseif ad==2
            spec.distr='T';
            spec.dparams=1;
            spec.dtheta0=10;
        else
            spec.distr='SkewT';
            spec.dparams=2;
            spec.dtheta0=[7;0];  
        end
        spec.vecsize = spec.mparams + spec.vparams + spec.dparams;
         if aa~=3
        ae=menu('define the PIT method','IFM','CML');
        if ae == 1
            spec.PIT = 'IFM';
        else
            spec.PIT = 'CML';
        end
        else
            spec.PIT = 'IFM';
        end
    end
    if aa~=1 
    if aa==2 || aa==3 %  possibly time varying copulas
            af = menu('define the copula family','t','Gaussian','Clayton','SJC');
            if af == 1
            spec.family = 't';
            elseif af == 2
            spec.family = 'Gaussian';
            elseif af == 3
            spec.family = 'Clayton';
            elseif af == 4
            spec.family = 'SJC';
            end
    elseif aa==4 || aa==5
            af = menu('define the copula family','t','Clayton','SJC');
            if af == 1
            spec.family = 't';
            elseif af == 2
            spec.family = 'Clayton';
            elseif af == 3
            spec.family = 'SJC';
            end
    end 
    end
    if aa~=1 && aa~=4 && aa~=5 % possibly time varying copulas
        if af == 1
            ag=menu('define the evolution of the copula parameter','static','DCC');
            if ag == 1
                spec.depspec = 'static';
                spec.ctheta0 = 10;
            elseif ag ==2
                spec.depspec = 'DCC';
                spec.ctheta0 = [12;.015;.98];
            end
        elseif af == 2
            spec.depspec = 'DCC';
            spec.ctheta0 = [.015;.98];
            display('for the Gaussian copula only the DCC specification is supported')
            display('the MLE of the correlation ofm the Gaussian copula is the sample')
            display('correlation of the transformed data')
        elseif af == 3 || af == 4
            ah=menu('define the evolution of the copula parameter','static','time varying');
            if ah == 1
                spec.depspec = 'static';
            elseif ah ==2
                spec.depspec = 'Patton';
            end
            if af == 3 && ah == 1 %static Clayton
                xxx = corr(data,'type','Kendall');
                spec.ctheta0 = xxx(1,2);
            elseif af == 3 && ah == 2 % tv Clayton
                xxx = corr(data,'type','Kendall');
                ctheta0 = xxx(1,2);
                spec.ctheta0 = [log(ctheta0/(1-ctheta0));0;0];
            elseif af ==4 && ah == 1 % static SJC
                spec.ctheta0 = .25*ones(2,1);
            elseif af == 4 && ah == 2 % tv SJC
                spec.ctheta0 = repmat([-1;-.5;.5],[2,1]);
            end
        end
        if aa == 3
        spec.theta0 = [repmat([spec.mtheta0;spec.vtheta0;spec.dtheta0],[n,1]);spec.ctheta0];
        end
    end
    if aa == 4
        spec.depspec = 'static';
        ai = menu('define copula vine decomposition','C - Vine', 'D - Vine');
        if ai == 1
            spec.decomposition = 'CVine';
        elseif ai == 2
            spec.decomposition = 'DVine';
        end
        if af == 1
            spec.ctheta0 = 10*ones(.5*n*(n-1),1);
        elseif af == 2
            spec.ctheta0 = .25*ones(.5*n*(n-1),1);
        elseif af == 3
            spec.ctheta0 = .25*ones(n*(n-1),1);
        end
    end
    if aa == 5
        spec.depspec = 'static';
        
        ai = menu('define copula vine decomposition','C - Vine', 'D - Vine');
        if ai == 1
            spec.decomposition = 'CVine';
        elseif ai == 2
            spec.decomposition = 'DVine';
        end
        if af == 1
            spec.ctheta0 = 10;
        elseif af == 2
            spec.ctheta0 = .25;
        elseif af == 3
            spec.ctheta0 = [.25;.25];
        end
    end
        
end
            
            
            
        