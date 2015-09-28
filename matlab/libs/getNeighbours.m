function [undirected, incoming, outgoing] = getNeighbours( dag, node )
%GETUNDIRECTEDEDGES Summary of this function goes here
%   return neighbours of a node where the direction is unknow.
    undirected = [];
    incoming = [];
    outgoing = [];
    for i = 1:length(dag(:, node))
        if dag(i, node) == 1 && dag(node, i) == 1
            % edge of (i, node) is undirected
            undirected = [undirected, i];
        end
        if dag(i, node) == 1 && dag(node, i) == 0
            % edge of (i, node) is incoming to node
            incoming = [incoming, i];
        end
        if dag(i, node) == 0 && dag(node, i) == 1
            % edge of (i, node) is outgoing from node
            outgoing = [outgoing, i];
        end
    end

end

