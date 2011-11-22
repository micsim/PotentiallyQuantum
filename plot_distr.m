function plot_distr(eigenvectors, distribution, v_vec)
    n = size(distribution,1);
    assert(size(eigenvectors,2) == n,...
        'The distribution vector must have the same size as the eigenvectors');

    x_values = (1:size(eigenvectors,1))';

    plot(x_values, eigenvectors * distribution, x_values, v_vec);
end
