    
    function E = a2e(A)
    % Convert Adjacency Matrix to Edge List
        [I,J] = find(A); %Find indices of nonzero elements.
        E = [I J];
    end