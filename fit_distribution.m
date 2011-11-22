function v = fit_distribution(eigenvectors, distribution)
    n = size(eigenvectors,2);
    m = size(distribution, 1);
    assert(size(eigenvectors,1) == m,...
        'The distribution vector must have the same size as the eigenvectors');

    v = zeros(n,1);

    for i = 1:n
        v(i) = eigenvectors(:,i)' * distribution;
    end
end
