% Tests a timedependent potential that pushes the particle.

n = 300;
L = 1e-5;
u = (1:n)';

n_t = 1000;
dt = 3e-12;

%d_orig = create_gauss_distr(n, L, 0.1, 0.03, 0 * pi);

v_max = 3.5e-5;
v_mat = zeros(n, n_t);
for i = 1:n_t
    min_pos = 1;
    max_pos = floor(i/n_t * (n-1)) + 1;
    v_mat(min_pos:max_pos,i) = v_max;

    % = max(min(floor((i/n_t + 0.5)*(n-1)) + 1, n), 0);
end

%psi = cn_evolution(v_mat, d_orig, L, dt, @get_kinetic_hamiltonian, 0.1);

%plot_timedep_slider(psi, n_t, dt, 1e-4*0.001*v_mat/(v_max^2 * norm(d_orig,inf)^2));


simulation = StepwiseSimulation(v_mat, 1, 0.1, 0.03, L, 0, 0.1);
simulation.dt = dt;
simulation.plot = true;