% v_func is the potential given as a function of x
% n is the number of points used for discretising the position.
% factor is h ("h-bar") K / 2m
% where K is the length of ...
function [eigenvalues, eigenvectors] = get_hamiltonian_eigenvectors(v_func, n, factor)
    assert(mod(n,2) == 0, 'n must be even!');

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
    H = diag(v_func(1:n)) + T;

    disp('The hamiltonian:');
    disp(H);

    % Find eigenvalues and eigenvectors:
    [eigenvectors,D] = eig(H);

    % The eigenvalues as a vector:
    eigenvalues = diag(D);
    % Vector of the values on the diagonal.
end

