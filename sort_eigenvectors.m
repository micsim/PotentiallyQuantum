function [eigenvalues, eigenvectors] = sort_eigenvectors(eigenvalues, V)
    n = size(eigenvalues,1);

    % Sort eigenvalues and store order:
    [eigenvalues, order] = sort(eigenvalues);

    % Create eigenvector array in advance:
    eigenvectors = zeros(size(V));

    % Also order eigenvectors:
    for j = 1:n
        eigenvectors(:,j) = V(:,order(j));
    end
end
