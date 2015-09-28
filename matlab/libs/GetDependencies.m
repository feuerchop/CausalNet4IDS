function dep = GetDependencies( graph, i, j )
%GETDEPENDENCIES Summary of this function goes here
%   return the dependencies set given an edge, corresponding to its markov
%   blanket of both nodes
        dep = [];
        seti = union(find(graph(i, :)), find(graph(:, i))');
        seti(find(seti == j)) = [];
        setj = union(find(graph(j, :)), find(graph(:, j))');
        setj(find(setj == i)) = [];
        dep = union(seti, setj);
end

