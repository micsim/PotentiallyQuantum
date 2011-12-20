function plot_timedep(psi, v_vec)
    n = size(psi,1);
    assert(n == size(v_vec,1));
    
    x_values = (1:n)' / n;

    plot(x_values, abs(psi).^2, x_values, v_vec);
end
