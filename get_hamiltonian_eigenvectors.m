% v_vec is the potential given as a vector of n values for each x.
% Where n is the number of points used for discretising the position.
% factor is h ("h-bar") K / 2m
% where K is the length of ...
function [eigenvalues, eigenvectors] = get_hamiltonian_eigenvectors(v_vec, factor)
    n = size(v_vec, 1);
    assert(size(v_vec,2) == 1, 'v_vec must be a column vector!');
    assert(mod(n,2) == 0, 'v_vec s size must be even!');

    % Build the hamiltonian:
    %  First the kinetic energy:
    %    The diagonal:
    T = -factor * 1/3 * (1 + 2/(n^2)) * eye(n);
    % An diagonal matrix of all the same values (that value times the identity).

    %    The other elements:
    for i = 1:n
        for j = 1:n
            if i ~= j
                T(i,j) = -factor * 2/(n^2) * (-1)^(j-1) / (sin(pi*(j-i)/n));
            end
        end
    end

    % The hamiltonian itself:
    H = diag(v_vec) + T;

    %disp('The hamiltonian:');
    %disp(H);

    % Find eigenvalues and eigenvectors:
    [V,D] = eig(full(H));
    
    % The eigenvalues as a vector:
    eigenvalues = diag(D);
    % Vector of the values on the diagonal.
    
    % Sort eigenvalues and store order:
    [eigenvalues,order] = sort(eigenvalues);
    
    % Create eigenvector array in advance:
    eigenvectors = zeros(n);
    
    % Also order eigenvectors:
    for j = 1:n
        eigenvectors(:,j) = V(:,order(j));
    end
end

