function candidates = findColliderCandidates( graph, p, q )
%FINDCOLLIDERCANDIDATES Summary of this function goes here
%   Find out collider candidates give an edge (p, q)
    candidates = [];
    if graph(p, q) ~= 1
        candidates = -1;
        return;
    else
        for i = 1:size(graph, 2)
            if graph(p, i) == 1 && graph(q, i) == 1
                candidates = [candidates; i];
            end
        end
    end

end

