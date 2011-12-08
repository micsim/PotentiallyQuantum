% v_vec is the potential given as a vector of n values for each x in eV.
% L is the length of the periodic box in position space in meters.
% mass is the mass of the particle relative to the electron mass.
% 
% eigenvalues is a sorted vector of the eigenvalues in eV.
% eigenvectors is a matrix whose columns are the eigenvectors in the same order
% as the corresponding eigenvalues.
function [eigenvalues, eigenvectors] = get_hamiltonian_eigenvectors(v_vec, L, eig_function, mass)
    n = size(v_vec, 1);
    % n is the number of points used for discretising the position.

    assert(size(v_vec,2) == 1, 'v_vec must be a column vector!');
    assert(mod(n,2) == 0, 'v_vec s size must be even!');

    if nargin < 4
        mass = 1;

        if nargin < 3
            eig_function = @(H) eig(full(H));
        end
    end

    h_bar = 1.054571506e-34; % [J*s]
    joule_to_eV = 6.24150974e18; % [eV/J]
    electron_mass = 9.1093829140e-31; % [kg]

    % factor is   h^2 ("h-bar") K^2  / 2m
    %           = h^2 ("h-bar") pi^2 / 2m(dx^2)
    %           = h^2 ("h-bar") pi^2 (n-1)^2 / (2m L^2)
    factor = joule_to_eV * h_bar^2 * pi^2 * (n-1)^2 / (2*mass*electron_mass * L^2); % [eV]

    % Build the hamiltonian:
    %  First the kinetic energy:
    %    The diagonal:
    T = factor * 1/3 * (1 + 2/(n^2)) * eye(n);
    % An diagonal matrix of all the same values (that value times the identity).

    %    The other elements:
    for i = 1:n
        for j = 1:n
            if i ~= j
                T(i,j) = factor * 2/(n^2) * (-1)^(j-i) / (sin(pi*(j-i)/n))^2;
            end
        end
    end

    % The hamiltonian itself:
    H = diag(v_vec) + T;

    %disp('The hamiltonian:');
    %disp(H);

    % Find eigenvalues and eigenvectors:
    [V,D] = eig_function(H);

    % normalize:
    V = sqrt(n/L)*V;
    
    % The eigenvalues as a vector:
    eigenvalues = diag(D);
    % Vector of the values on the diagonal.

    % Sort eigenvalues and store order:
    [eigenvalues, eigenvectors] = sort_eigenvectors(eigenvalues, V);
    % TODO Optimize: Not needed for using eigs as eig_function!
end
