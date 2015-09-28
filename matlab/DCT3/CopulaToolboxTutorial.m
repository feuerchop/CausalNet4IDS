function CopulaToolboxTutorial
clc
display('************** TOOLBOX TUTORIAL **********************')
display(' ')
display('For more information about the toolbox read the attached pdf file (ref1).');
display('This is just a quick tutorial.')
display(' ')
display('author: Manthos Vogiatzoglou, University of Macedonia, Thessaloniki Greece')
display(' ')
display('contact at: vogia@yahoo.com')
display(' ')
display(' --------- press enter to continue -----------------')
display(' ')
pause
clc
aa = menu('what do you want to learn about?','troubleshooting','how to use the toolbox');
if aa==1
    display('If during the use of the toolbox the procedure crashes, send your data to the')
    display('to the author, at: vogia@yahoo.com. Before doing so look at the following.')
    display('In the Clayton or SJC vines the parameters, Kendalls tau for the former and')
    dislay(' tail dependence for the latter belong in theory to (0,1)')
    display('however for loosely dependent data set, even moderately large values of')
    display('the parameters may cause the log - likelihood to go to -Inf')
    display('in that case, change the upper bound of the fmincon function to something')
    display('smaller and run the optimization routine again')
    display('when the standard errors are inf or NaN, some of the parameters is at the bound')
    display('concider changing the bounds and try again. The bounds are created within')
    display('the function CreateFminconConstraints. If the problem appears with the fminunc')
    display('it is probably because a+b>1 in a GARCH or DCC model. Try changing starting')
    display('values or experiment with the function RescaleParameters that implements')
    display('the reparametrization of the parameters.')
    display(' ')
    display('press enter to continue')
    pause
elseif aa==2
    bb=menu('What type of data do you have?','unfiltered (e.g stock returns)','series with iid margins','matrix with U(0,1) margins');
    if bb==1
        display('to estimate a copula, your data should have uniform U(0,1) margins.')
        display('first you have to decide in how many steps you want to fit your model')
        display('the one step method is more efficient but more computationaly involved')
        display('two step methods are computationaly efficient')
        display('press enter to continue.')
        pause
        display('For the 2 - step procedure run modelspec and in "what do you want to estimate?"')
        display('select choices one: GARCH model for each series')
        display('Choice 3 fits a copula model in one step')
        display('If you choose 1, your margins are assumed as GARCH models')
        display('PIT method defines how the iid residuals are transformed to uniform')
        display('If CML is selected the transformation is being made semiparametricaly')
        display('with the empiricalCDF function, else the distribution assumed for the GARCH model is used')
        display('After the modelspec, call the fitModel function to fit the defined model to the series')
        display('The output "udata" contains your iid uniform series.')
        display('Now you can fit a copula to udata')
        display(' ')
    elseif bb==2
        display('to estimate a copula, your data should have uniform U(0,1) margins.')
        display('If you know the distribution Fi, of the i - margin, you can transform')
        display('it to uniform by applying Fi to the i - margin:')
        display(' ')
        display('           if x~F then F(x)~U(0,1)               ')
        display(' ')
        display(' --------- press enter to continue -----------------')
        display(' ')  
        pause
    elseif bb==3
        display('*********** Parameter estimation ********************** ')
        display(' ')
        display('the toolbox is written to estimate via maximum likelihood the model parameters')
        display('run modelspec to define you copula or copula vine and then call')
        display('fitModel to estimate the parameters. Do not forget that the input data')
        display('should consist of uniform iid margins!!!')
        display('Consider also that for the one step copulas or copula vines the use')
        display('of good starting values is essential. For the one step models, estimate')
        display('the model in two steps first and then use these estimates for starting values.')
        display('For the vines, first choose to estimate "copula vine sequentially"')
        display('and use the output as starting values for the copula vine')
        display(' ')
    end
    
    bc=menu('display F.A.Q. menu?','yes','no');
    if bc==1
    cc=menu('FAQ','which specifications are supported by the copula?','what solver should I choose, fmincon or fminunc?');
    if cc==1
        display('******* Supported Copula specifications *************')
        display(' ')
        display('the toolbox support the following copulas:')
        display('the Clayton Copula, the SJC copula, the t - Copula and the Gaussian Copula,')
        display('the first two are supported only for the bivariate case. For the latter two')
        display('there is no dimensional constraint. In all cases the copula parameters')
        display('can be chosen by the user to be static or time varying.')
        display(' ')
        display(' --------- press enter to continue -----------------')
        display(' ')
        pause
        display('for time varying copulas the following specifications are supported:')
        display('In the t and Gaussian copula, the dependence parameter (correlation) can')
        display('follow dynamics similar to the DCC model of Engle.')
        display('For the Clayton and SJC copulas, the dependence parameter')
        display('follows a specification introduced in Patton (ref:2).')
        display(' ')
        display(' --------- press enter to continue -----------------')
        display(' ')
        pause
        display('In the static t copula, the correlation matrix equals the sample correlation')
        display('thus only the degree of freedom is estimated. The MLE of the Gaussian copula ')
        display('parameter is the sample correlation therefore there is no need to estimate')
        display('a static Gaussian Copula.')
        display(' ')
        display(' --------- press enter to continue -----------------')
        display(' ')
        pause
        display('for copula vines, both D - Vines and Canonical vines are supported')
        display('with all bivariate copulas in the cascade to belong to one of the following families:')
        display('t copulas, Clayton copulas, or SJC copulas')
    elseif cc==2
        display('******************** Solver choice **********************')
        display(' ')
        display('the problems are constrained by nature therefore fmincon is the obvious choice.')
        display('however fmincon does not provide a good approximation of the Hessian matrix.')
        display('therefore if one wants to calculate numerical standard error, fminunc should be used.')
    end
    end
    dd=menu('display references?','yes','no');
    display(' ')
    if dd==1
        display('ref1:Financial Copula toolbox. (submitted to the journal of statistical software)')
        display(' ')
        display('ref2:Modelling asymmetric exchange rate dependence. International economic review, 47 - 2, 2006')
        display(' ')
    end
end


