function plot_distr(eigenvectors, distribution, v_vec)
    n = size(eigenvectors,1);
    assert(size(distribution,2) == 1, 'distribution must be a column vector.');
    assert(size(v_vec,2) == 1, 'v_vec must be a column vector.');
    assert(n == size(v_vec,1),...
        'The potential vector must have the same size as the eigenvectors.');
    assert(size(eigenvectors,2) == size(distribution,1),...
        'The distribution vector must have the same size as the number eigenvectors.');

    x_values = (1:n)' / n;
    y_values = eigenvectors * distribution;

    plot(x_values, abs(y_values).^2, x_values, real(y_values).^2,...
         x_values, imag(y_values).^2, x_values, v_vec);
    %plot(x_values, abs(y_values).^2, x_values, v_vec);
end
