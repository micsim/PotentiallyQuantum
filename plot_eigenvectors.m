function plot_eigenvectors(eigenvectors, num)
    if nargin < 2
        num = inf;
    end

    num = min(num,size(eigenvectors,2));

    eigenvectors = eigenvectors(:,1:num);

    probabilities = arrayfun(@(ele) ele*ele', eigenvectors);

    x_values = (1:size(eigenvectors,1))';

    plot_arg = cell(num,1);
    for i = 1:num
        plot_arg{2*i - 1} = x_values;
        plot_arg{2*i}     = probabilities(:,i);
    end

    plot(plot_arg{:});
end

