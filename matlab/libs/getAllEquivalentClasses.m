function equivDags= getAllEquivalentClasses( pdag, equivDags, undirectedEdges )
%GETALLEQUIVALENTCLASSES Summary of this function goes here
%   from a partially directed acyclic graph to obtain all the equivalent
%   classes in a recursive method
    if ~isempty(undirectedEdges)
        pdag1 = pdag;
        pdag2 = pdag;
        pdag1(undirectedEdges(1, 2), undirectedEdges(1, 1)) = 0;
        pdag2(undirectedEdges(1, 1), undirectedEdges(1, 2)) = 0;
        undirectedEdges1 = undirectedEdges;
        undirectedEdges2 = undirectedEdges;
        undirectedEdges1(1, :) = [];
        undirectedEdges2(1, :) = [];
        
        [pdag1, undirectedEdges1] = ConstrainPropagation(pdag1, undirectedEdges1);
        [pdag2, undirectedEdges2] = ConstrainPropagation(pdag2, undirectedEdges2);
        if isempty(undirectedEdges1)
            equivDags{length(equivDags) + 1, 1} = pdag1;
        else
            equivDags = getAllEquivalentClasses(pdag1, equivDags, undirectedEdges1);
        end
        if isempty(undirectedEdges2)
            equivDags{length(equivDags) + 1, 1} = pdag2;
        else
            equivDags = getAllEquivalentClasses(pdag2, equivDags, undirectedEdges2);
        end
    else
        equivDags{length(equivDags) + 1, 1} = pdag;
    end
end

