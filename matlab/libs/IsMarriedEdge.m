function [ flag, colliders ] = IsMarriedEdge( graph, corr, p, q, alpha )
%ISMARRIEDEDGE Summary of this function goes here
%   Find out if an edge is a married edge brought by colliders, p and q are
%   the nodes of the edge
    flag = -1;
    colliderCandidates = findColliderCandidates(graph, p, q);
    if colliderCandidates == -1
        return
    end
    dependencies = GetDependencies(graph, p, q);
    for i = 1:length(colliderCandidates)
        % maximal so many colliders for one edge
        % loop breaks until find real colliders
        % find combination of colliders of size i
        colliderComb = nchoosek(colliderCandidates, i);
        for j = 1:size(colliderComb, 1)
            % traverse all the collider combinations
            colliders = colliderComb(j, :);
            setWithoutCollider = dependencies;
            for k = 1:length(colliders)
                setWithoutCollider( find(setWithoutCollider == colliders(k)) ) = [];
            end
            setWithoutCollider = [p, q, setWithoutCollider];
            setWithCollider = [p, q, dependencies];
%             part_invcorr = ScaleMatrix(inv( corr(setWithoutCollider, setWithoutCollider) )) ;
            part_invcorr = ScaleMatrix((corr(setWithoutCollider, setWithoutCollider) + 1e-5*eye(length(setWithoutCollider)))\eye(length(setWithoutCollider)));
            %part_invcorr = inv( corr(setWithoutCollider, setWithoutCollider) );
            if abs(part_invcorr(1, 2)) < alpha
            %if ~FisherSignificanceTest(N, length(corr)-2, alpha, part_invcorr(1, 2))
                full_invcorr = ScaleMatrix((corr(setWithCollider, setWithCollider)+1e-5*eye(length(setWithCollider)))\eye(length(setWithCollider)));
                %full_invcorr = inv( corr(setWithCollider, setWithCollider) );
                if abs(full_invcorr(1, 2)) > alpha
                %if FisherSignificanceTest(N, length(corr)-2, alpha, full_invcorr(1, 2))
                    flag = 1;
                    return;
                end
            end
        end
    end
    if flag == -1
        % no colliders are found to be true
        colliders = [];
    end
end

