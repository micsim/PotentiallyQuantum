function [eigenvectors] = evolve_eigenvectors(eigenvalues, eigenvectors, t)
    h_bar = 1;
    factor = -i/h_bar * t;
    for j = 1:size(eigenvalues,1)
        eigenvectors(:,j) = exp(factor * eigenvalues(j)) * eigenvectors(:,j);
    end
end

