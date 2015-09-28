function  result = sampleBN( bnet, count )
%GENERATEDATAFROMBN Summary of this function goes here
%   Detailed explanation goes here
    i = 1;
    result = [];
    while i <= count 
        sample = sample_bnet(bnet);
        for j = 1:length(sample)
            result(i,  j) = sample{j};            
        end
        i = i + 1;      
    end
end

