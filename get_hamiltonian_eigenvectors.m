% v_vec is the potential given as a vector of n values for each x in eV.
% L is the length of the periodic box in position space in meters.
% mass is the mass of the particle relative to the electron mass.
% 
% eigenvalues is a sorted vector of the eigenvalues in eV.
% eigenvectors is a matrix whose columns are the eigenvectors in the same order
% as the corresponding eigenvalues.
function [eigenvalues, eigenvectors,H] = get_hamiltonian_eigenvectors(v_vec, L, eig_function, kin_h_function, varargin)
    n = size(v_vec, 1);
    % n is the number of points used for discretising the position.

    if nargin < 4
        kin_h_function = @get_kinetic_hamiltonian;

        if nargin < 3
            eig_function = @(H) eig(full(H));
        end
    end

    assert(size(v_vec,2) == 1, 'v_vec must be a column vector!');

    % The hamiltonian itself:
    H = diag(v_vec) + kin_h_function(n,L,varargin{:});

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
