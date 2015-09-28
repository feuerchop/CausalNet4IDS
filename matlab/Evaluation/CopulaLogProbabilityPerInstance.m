function lp = CopulaLogProbabilityPerInstance(trainData, testData, ctype, rho)
    lp = 0;
    if strcmp(ctype, 'gauss')
        for i = 1:size(testData, 1)
            pr = gaussian_copula_bn_jpdf(testData(i, :), trainData, rho);
            if pr <=0 
                pr = 0.00001;
            end
            fprintf('jpdf: %f\n', pr);
            lp = lp + log(pr);
        end
        lp = lp/size(testData, 1);
        fprintf('Avg Logprob on testdata: %f', lp);
    else
        lp = -1;
        display('Copula Type Not Suported Yet.');
    end
end