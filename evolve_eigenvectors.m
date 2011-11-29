% t is the time in seconds
function [eigenvectors] = evolve_eigenvectors(eigenvalues, eigenvectors, t)
    h_bar = 1.054571506e-34; % [J*s]
    joule_to_eV = 6.24150974e18; % [eV/J]

    factor = -i/(h_bar * joule_to_eV) * t;
    for j = 1:size(eigenvalues,1)
        eigenvectors(:,j) = exp(factor * eigenvalues(j)) * eigenvectors(:,j);
    end
end
