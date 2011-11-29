% v_vec is the potential given as a vector of n values for each x.
% Where n is the number of points used for discretising the position.
% factor is h^2 ("h-bar")/2m(dx^2)
% where dx is the step-width of position-space.
% This function uses the difference quotient to calculate solutions
% of the one-dimensional time-independent Schrödinger equation.
function [eigenvalues, eigenvectors] = get_hamiltonian_eigenvectors_dq(v_vec, factor)
    n = size(v_vec, 1);
    assert(size(v_vec,2) == 1, 'v_vec must be a column vector!');

    % Build the hamiltonian:
    %  First the kinetic energy:
    %    The diagonal:
    T = factor * 2* eye(n);
    % An diagonal matrix of all the same values (that value times the identity).

    %    The other elements:
    for i = 2:n
        T(i,i-1) = -factor;
        T(i-1,i) = -factor;
    end
    
    % Add border values for periodic boundary conditions
    T(1,n) = -factor;
    T(n,1) = -factor;

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

