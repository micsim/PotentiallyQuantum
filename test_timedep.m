joule_to_eV = 6.24150974e18; % [eV/J]
h_bar = 1.054571506e-34; % [J*s]
h_bar = h_bar * joule_to_eV; % covert to [eV*s]

n = 300;
u = (1:n)';
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 3.5e-5;

n_t = 200;
dt = 1e-10;

psi = zeros(n, n_t);

k = 1e7 * pi;

d_orig = sqrt(normpdf(u, 0.1*n, 0.03*n));
for j = 1:n
    pos = j/n * L;
    d_orig(j) = exp(1i*k*pos) * d_orig(j);
end
psi(:,1) = d_orig/norm(d_orig);

T = get_kinetic_hamiltonian(n,L);

for i = 1:n_t-1
    %H = diag(v_vec(:,i)) + T;
    H = diag(v) + T;
    H2 = 1i/h_bar * H * dt / 2;

    psi(:, i+1) = ((eye(n) + H2)\(eye(n) - H2)) * psi(:,i);
    % DOES NOT WORK: TODO why?
    %tmp = (eye(n) + 2*H2) * psi(:,i);
    %psi(:, i+1) = tmp/norm(tmp);
end

plot_timedep_slider(psi, n_t, dt, 1e1*0.001*v/(norm(v,inf) * norm(d_orig,inf)^2));
