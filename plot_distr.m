function plot_distr(eigenvectors, distribution, v_vec)
    assert(size(eigenvectors,2) == size(distribution,1),...
        'The distribution vector must have the same size as the eigenvectors');

    x_values = (1:size(eigenvectors,1))';
    y_values = eigenvectors * distribution;

    plot(x_values, abs(y_values).^2, x_values, real(y_values),...
         x_values, imag(y_values), x_values, v_vec);
end
