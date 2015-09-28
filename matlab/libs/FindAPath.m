function [flag, path] = FindAPath( pdag, p, q )
%FINDAPATH Summary of this function goes here
%   dijisktra algorithm
    flag = 0; path = [];
    dist = ones(1, size(pdag, 1));
    previous = ones(1, size(pdag, 1));
    nodes = [1:size(pdag, 1)];
    dist(:) = inf;
    previous(:) = -1;
    dist(p) = 0;
    while ~isempty(nodes)
        n = find( dist(nodes) == min(dist(nodes)) );
        n = nodes(n(1));
        if dist(n) == inf
            break;
        else
            nodes( find(nodes == n) ) = [];
            [undirected, incoming, outgoing] = getNeighbours( pdag, n );
            for i = 1:length(outgoing)
                v = outgoing(i);
                if (dist(n) + 1) < dist(v)
                    dist(v) = dist(n) + 1;
                    previous(v) = n;
                end
            end
        end
    end
    if dist(q) ~= inf
        flag = 1;
        path = [q];
        while previous(q) ~= -1
            path = [path, previous(q)];
            q = previous(q);
        end
    end
    
    
end

