
    function D = ve2d(V,E)
    % Compute Euclidean Distance for Edges
        VI = V(E(:,1),:);
        VJ = V(E(:,2),:);
        D = sqrt(sum((VI - VJ).^2,2));
    end