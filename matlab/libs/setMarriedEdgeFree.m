function dag = setMarriedEdgeFree( dag,  i, j, colliders)
%SETMARRIEDEDGEFREE Summary of this function goes here
%   suppose the edge (i, j) is a married edge brought by colliders, then
%   set it free and configure the correct direction to colliders
    dag(i, j) = 0;
    dag(j, i) = 0;
    for k = 1:length(colliders)
        collider = colliders(k);
        dag(i, collider) = 1;
        dag(collider, i) = 0;
        dag(j, collider) = 1;
        dag(collider, j) = 0;
    end

end

