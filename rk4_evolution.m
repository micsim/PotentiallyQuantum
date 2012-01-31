%Runge-Kutta-Implementation: Achtung: dt = Zeitschritt der
%Potentialdiskretisierung, psi-diskretisierung = 2 * dt!!
function psi = rk4_evolution(v_mat, psi_0, L, dt, kin_h_function, mass)
    joule_to_eV = 6.24150974e18; % [eV/J]
    h_bar = 1.054571506e-34; % [J*s]
    h_bar = h_bar * joule_to_eV; % covert to [eV*s]

    n = size(v_mat, 1);
    % n is the number of points used for discretising the position.

    n_t = size(v_mat, 2)/2;
    % n_t is the number of timesteps. Runge-Kutta: V has to be evaluated
    % between two evaluations of psi

    assert(size(psi_0, 1) == n, 'psi_0 must have the same dimension as the potential in one point in time.');
    assert(size(psi_0, 2) == 1, 'psi_0 must be a column vector!');

    if nargin < 6
        mass = 1;
        if nargin < 5
            kin_h_function = @get_kinetic_hamiltonian;
        end
    end


    psi = zeros(n, n_t);
    psi(:,1) = psi_0;

    T = kin_h_function(n,L,mass);
    
    h = dt*2 %length of one timestep in the psi discretization
     
    for j = 1:n_t-1
        
        psi_j = psi(:,j);
        H1 = diag(v_mat(:,2*j-1)) + T;
        k1 = 1/(1i*h_bar)*H1*psi_j;
        %display(size(k1,1));
        %display(k1);
        
        H2 = diag(v_mat(:,2*j)) + T;
        k2 = 1/(1i*h_bar)*H2*(psi_j+h/2*k1);
        k3 = 1/(1i*h_bar)*H2*(psi_j+h/2*k2);
        
        H3 = diag(v_mat(:,2*j+1)) + T;
        k4 = 1/(1i*h_bar)*H3*(psi_j+h*k3);
        
        psi(:,j+1) = psi_j+h/6*(k1 + 2*(k2+k3) + k4);
        
        %if j==4
         %   break;
        %end
        
        %H2 = 1i/h_bar * H * dt / 2;

        %psi(:, i+1) = ((eye(n) + H2)\(eye(n) - H2)) * psi(:,i);
        % DOES NOT WORK: TODO why?
        %tmp = (eye(n) + 2*H2) * psi(:,i);
        %psi(:, i+1) = tmp/norm(tmp);
    end
end
