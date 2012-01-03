function psi = cn_evolution(v_mat, psi_0, L, dt, kin_h_function, varargin)
    joule_to_eV = 6.24150974e18; % [eV/J]
    h_bar = 1.054571506e-34; % [J*s]
    h_bar = h_bar * joule_to_eV; % covert to [eV*s]

    n = size(v_mat, 1);
    % n is the number of points used for discretising the position.

    n_t = size(v_mat, 2);
    % n_t is the number of timesteps.

    assert(size(psi_0, 1) == n, 'psi_0 must have the same dimension as the potential in one point in time.');
    assert(size(psi_0, 2) == 1, 'psi_0 must be a column vector!');

    if nargin < 5
        kin_h_function = @get_kinetic_hamiltonian;
    end


    psi = zeros(n, n_t);
    psi(:,1) = psi_0;

    T = kin_h_function(n,L,varargin{:});

    for i = 1:n_t-1
        H = diag(v_mat(:,i)) + T;
        H2 = 1i/h_bar * H * dt / 2;
    
        psi(:, i+1) = ((eye(n) + H2)\(eye(n) - H2)) * psi(:,i);
        % DOES NOT WORK: TODO why?
        %tmp = (eye(n) + 2*H2) * psi(:,i);
        %psi(:, i+1) = tmp/norm(tmp);
    end
end
