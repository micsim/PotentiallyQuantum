function plot_amplitudes(eigenvectors, eigenvalues, num)
    assert(size(eigenvalues,2) == 1, 'eigenvalues must be a column vector.');
    assert(size(eigenvalues,1) == size(eigenvectors,1),...
        'The number of eigenvalues must equal the size of the eigenvectors.');

    if nargin < 3
        num = inf;
    end

    num = min(num,size(eigenvectors,2));

    eigenvectors = eigenvectors(:,1:num);

    probabilities = arrayfun(@(ele) abs(ele)^2, eigenvectors);

    x_values = (1:size(eigenvectors,1))';

    plot_arg = cell(num,1);
    for i = 1:num
        plot_arg{2*i - 1} = x_values;
        plot_arg{2*i}     = real(eigenvectors(:,i));
    end

    plot(plot_arg{:});
end

