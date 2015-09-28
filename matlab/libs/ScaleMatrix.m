function s = ScaleMatrix( m )
%SCALEMATRIX Summary of this function goes here
%   Detailed explanation goes here
    for i = 1:size(m, 1)
        for j = 1:size(m, 2)
            if i == j
                s(i, j) = 1;
            else
                s(i, j) = m(i, j)/sqrt(m(i, i)*m(j, j));
            end
        end
    end

end

