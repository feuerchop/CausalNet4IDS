function flag = checkDiamondStructure( dag, p, q )
%CHECKDIAMONDSTRUCTURE Summary of this function goes here
% check if there is a diamond structure between p, q which means
% there are two paths from p to q, p--w-->y and q--z-->y, p and q are
% connected as well. Note that p--w and p--z are not directed. w and z are
% not connected. In this case, p should point to q unless it will introduce
% new V-structure.
    flag = 1;
    [undirected_p, incoming, outgoing] = getNeighbours(dag, p);
    [undirected, incoming_q, outgoing] = getNeighbours(dag, q);
    vertexset = [];
    if length(undirected_p) < 3 || isempty(find(undirected_p == q))
        flag = 0;
        return;
    else
        for i = 1:length(undirected_p)
            if ~isempty( find( incoming_q == undirected_p(i) ) )
                vertexset = [vertexset, undirected_p(i)];
            end
        end
        if length(vertexset) < 2
            flag = 0;
            return;
        else
            for i = 1:length(vertexset)-1
                for j = i+1:length(vertexset)
                    if dag(vertexset(i), vertexset(j)) == 0 && dag(vertexset(j), vertexset(i)) == 0 
                        % w and z are not connected, this is a diamond
                        % structure, p should point to q
                        flag = 1;
                        return;
                    end
                end
            end
        end
    end
        
end

