function [cpdag, undirectedEdges]= ConstrainPropagation( pdag, undirectedEdges )
%CONSTAINPROPAGATION Summary of this function goes here
%   Propagate constains to obtain a complete partially acyclic graph
%   The order of rules is very critical. Acyclicity should be always firstly
%   considered unless there might be a cycle in final result
    % a flag indicating wether the graph is updated
    changed = 1;  directable = 0;
    while changed
        % if graph has any changes, then it should not to stop propagation
        changed = 0; directable = 0;
        % preserve acyclicity
        pointer = 1;
        LenOfUndirectedEdges = size(undirectedEdges, 1);
        while pointer <= LenOfUndirectedEdges           
            i = undirectedEdges(pointer, 1);
            j = undirectedEdges(pointer, 2);
            % path from i to j
            %reachableList = -1*ones(1, size(pdag, 1));
            if FindAPath(pdag, i, j) > 0
            % the direct i to j, unless there would be a cycle
                pdag(j, i) = 0;
                changed = 1;
                directable = 1;
                undirectedEdges(pointer, :) = [];
                LenOfUndirectedEdges = size(undirectedEdges, 1);
                break;
            elseif FindAPath(pdag, j, i) > 0
                pdag(i, j) = 0;
                changed = 1;
                directable = 1;
                undirectedEdges(pointer, :) = [];
                LenOfUndirectedEdges = size(undirectedEdges, 1);
                break;
            else
                pointer = pointer + 1;
            end          
        end   
        if directable 
            continue;
        end
                    
        % No new V-structure
        pointer = 1; 
        LenOfUndirectedEdges = size(undirectedEdges, 1);
        while pointer <= LenOfUndirectedEdges 
             i = undirectedEdges(pointer, 1);
             j = undirectedEdges(pointer, 2);
             directable = 0;
             [undirected, incoming, outgoing] = getNeighbours( pdag, i );
             if ~isempty(incoming)
                 for k = 1:length(incoming)
                     if pdag(incoming(k), j) == 0 && pdag(j, incoming(k)) == 0
                         pdag(j, i) = 0;
                         directable = 1;
                         break;
                     end
                 end
             else
                 [undirected, incoming, outgoing] = getNeighbours( pdag, j );
                 if ~isempty(incoming)
                    for k = 1:length(incoming)
                         if pdag(incoming(k), i) == 0 && pdag(i, incoming(k)) == 0
                            pdag(i, j) = 0;
                            directable = 1;
                            break;
                         end
                    end
                 end
             end
             if directable
                changed = 1;
                undirectedEdges(pointer, :) = [];
                LenOfUndirectedEdges = size(undirectedEdges, 1);    
                break;
             else
                 pointer = pointer + 1;
             end
        end
        if directable
            continue;
        end
        
        % form three-fork V structure with married parents
        pointer = 1;
        LenOfUndirectedEdges = size(undirectedEdges, 1);
        while pointer <= LenOfUndirectedEdges           
            i = undirectedEdges(pointer, 1);
            j = undirectedEdges(pointer, 2);     
            directable = 0;
           % path from i to j
           if checkDiamondStructure(pdag, i, j)
           % it forms a diamond structure between i and j, so
           % point i to j.
                pdag(j, i) = 0;
                changed = 1;
                directable = 1;
                undirectedEdges(pointer, :) = [];
               LenOfUndirectedEdges = size(undirectedEdges, 1);
               break;
           elseif checkDiamondStructure(pdag, j, i)
               pdag(i, j) = 0;
               changed = 1;
               directable = 1;
               undirectedEdges(pointer, :) = [];
               LenOfUndirectedEdges = size(undirectedEdges, 1);    
               break;
           else
                pointer = pointer + 1;
           end    
        end
        if directable
            continue;
        end
    end
    cpdag = pdag;
end

