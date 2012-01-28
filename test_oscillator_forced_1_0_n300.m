n = 300;
u = (1:n)';
L = 1e-5;
n_t = 8000;

omega = 1e9; % [rad/s]
dt = 16*(2*pi/(omega*n_t));%3e-11;


electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

omega_pot = 1*omega;  %angular frequency of the harmonically oscillating potential
ampl_pot = 0.03;        %amplitude of the potential oscillation, relative to L



d_orig = sqrt(normpdf(u, (0.5)*n, 0.03*n));


%psi = cn_evolution(v_mat, d_orig, L, dt);
load('pot_1-0_res.mat');
load('psi_1-0_res_n300.mat');


plot_timedep_slider(psi, n_t, dt, 4e0*0.001*v_mat/(norm(v_mat(:,1),inf) * norm(d_orig,inf)^2));

exp_values = zeros(n_t,1);

for l = 1:n_t
    exp_values(l,1) = u'*(abs(psi(:,l)).^2)/(norm(psi(:,l)))^2;
end
plot(exp_values);
