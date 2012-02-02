n = 1000;
L = 1e-5;

u = (1:n)';
a = floor(0.250*n);
b = floor(0.50*n);

% Initialise:
v = zeros(n,1);
v_s = v;              % sawtooth
v_t = v;              % triangle
v_b = v;              % block
v_wa = v;             % wall
v_we = v; %#ok<NASGU> % well
v_ho = v; %#ok<NASGU> % harmonic oscillator


v_2 = (u-a);
for i = 1:b-a
    v_2(b+i) = v_2(b-i);
end

% sawtooth potential:
v_s_width = floor(n/10);
for i = 1:n
    v_s(i) = mod(i,v_s_width);
end

% triangle:
v_t_width = floor(n/10);
for i = 1:n
    if mod(i,v_t_width*2) < v_t_width
        v_t(i) = mod(i,v_t_width)/v_t_width;
    else
        v_t(i) = 1 - mod(i,v_t_width)/v_t_width;
    end
end

% block:
v_b_width = floor(n/10);
for i = 1:n
    if mod(i, v_b_width*2) < v_b_width
        v_b(i) = 0;
    else
        v_b(i) = 1;
    end
end

% wall:
v_wa(floor(0.25*n):floor(0.5*n)) = 1;

% well:
v_we = 1 - v_wa;

% harmonic oscillator:
omega = 5e8; % [rad/s]
electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

v_ho = (((u/n - 0.5).*L).^2).*(0.5*electron_mass*omega^2*joule_to_eV);



% Set v to desired potential:
v = v_b;
%v(a:b) = v_2(a:b);

% Rescale v:
v_factor = 1e-5;
%v_factor = 8e-6;
v = v_factor * v/max(v);



L = 1e-5;
k = 5e6 * pi;

sim = EVSimulation(v,  5e-9,  0.1,    0.03, L, 5e6 * pi);
%                  V, t_max, d_mu, d_sigma, L,        k

% Compute and plot:
%sim.plot = true;