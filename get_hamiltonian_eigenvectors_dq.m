% v_vec is the potential given as a vector of n values for each x in eV.
% L is the length of the periodic box in position space in meters.
% mass is the mass of the particle relative to the electron mass.
% This function uses the difference quotient to calculate solutions
% of the one-dimensional time-independent Schrödinger equation.
function [eigenvalues, eigenvectors] = get_hamiltonian_eigenvectors_dq(v_vec, L, mass)
    n = size(v_vec, 1);
    % n is the number of points used for discretising the position.

    assert(size(v_vec,2) == 1, 'v_vec must be a column vector!');

    if nargin < 3
        mass = 1;
    end

    h_bar = 1.054571506e-34; % [J*s]
    joule_to_eV = 6.24150974e18; % [eV/J]
    electron_mass = 9.1093829140e-31; % [kg]

    % factor is   h^2 ("h-bar")/ 2m(dx^2)
    %           = h^2 ("h-bar") (n-1)^2 / (2m L^2)
    % where dx = L/(n-1) is the step-width of position-space.
    factor = joule_to_eV * h_bar^2 * (n-1)^2 / (2*mass*electron_mass * L^2); % [eV]

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

    % normalize:
    V = sqrt(n/L)*V;
    
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

