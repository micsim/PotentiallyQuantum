function v = fit_distribution(eigenvectors, distribution, L)
    n = size(eigenvectors,1);
    m = size(eigenvectors,2);
    assert(size(distribution, 1) == n,...
        'The distribution vector must have the same size as the eigenvectors');

    v = zeros(m,1);

    for i = 1:m
        v(i) = eigenvectors(:,i)' * distribution * L/n;
    end
end
