% Create a gauss distribution
function distribution = create_gauss_distr(n, L, mu, sigma, k)
    distribution = sqrt(normpdf((1:n)', mu*n, sigma*n));
    if k ~= 0
        %for j = 1:n
        %    pos = j/n * L;
        %    distribution(j) = exp(1i*k*pos) * distribution(j);
        %end
        distribution(1:n) = exp(1i*k*L*((1:n)'/n)) .* distribution(1:n);
        % Does the same as the above.
    end
    distribution = distribution/norm(distribution);
end