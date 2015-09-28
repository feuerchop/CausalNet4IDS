function [dag, colliderList] = detriangulation( graph, corr, alpha )
%FINDTRIANGULAR Summary of this function goes here
%   find out hanged triangulars from given node
    dag = graph;
    colliderList = [];
    for p = 1:size(dag, 2)
        for q=p+1:size(dag, 1)            
                if dag(p, q) ~= 0 && dag(q, p) ~= 0
                    % each unknown edge should be checked
                    [flag, colliders] = IsMarriedEdge(dag, corr, p, q, alpha);
                    for c = 1:length(colliders)
                            % gurantee there is no circle, unless drop this
                            % V-structure
                        if FindAPath(dag, colliders(c), p) || FindAPath(dag, colliders(c), q)
                                flag = -1;
                                break;
                        end
                    end
                    if flag ~= -1
                        dag = setMarriedEdgeFree(dag, p, q, colliders);
                        for k = 1:length(colliders)
                            % put collider into a global list
                            if isempty( find(colliderList == colliders(k)) )
                                colliderList = [colliderList, colliders(k)];
                            end
                        end
                    end                 
                end
        end
    end
end

