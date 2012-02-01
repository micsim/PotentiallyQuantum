classdef Simulation < handle
    %SIMULATION Parent class for EVSimulation, StepwiseSimulation.
    %       Simulate a particle in 1D with periodic boundary conditions.
    
    properties
        V = [];
        % Potential timeindependent -> n x 1       Matrix
        %           timedependent   -> n x n_t - 1 Matrix
        L = 1e-5;
        % Length of simulated space.
        k = 0;
        % k to start with.
        d = [];
        % The probability distribution of the particle at the starting
        % time.
        d_mu = 0.5;
        d_sigma = 0.1;
        % If d_mu and d_sigma are not empty, compute d from them in
        % reevolve as a normal distribution.
        mass = 1;
        % The mass of the particle relative to the electron mass.
        kin_h_function = @get_kinetic_hamiltonian;
        
        plot = false;
        % Wether the class is mainting a plot at the moment.
        t_max;
        % Maximum time of the slider.
    end
    
    properties (SetAccess = protected)
        T = [];
        % Kinetic hamiltonian.
    end
    
    properties (Dependent = true, SetAccess = protected)
        n;
        % Number of points the space is discretised in.
    end
    
    properties (Access = protected)
        recompute_pending    = false;
        redistribute_pending = false;
                
        % == GUI ==
        V_factor;
        % Factor for displaying the potential with a reasonable scale in the
        % y direction.
        
        x_plot;
        % Window for plotting the particle in position space.
        k_plot;
        % Window for plotting the particle in k space.
        slider;
        % Slider in the x_plot window for changing the observed time.
    end
    
    methods
        function o = Simulation(V, t_max, d_mu, d_sigma, L, k, mass, kin_h_function)
            o.V = V;
            o.t_max = t_max;
            
            if nargin > 2
                o.d_sigma = d_sigma;
                o.d_mu = d_mu;
                
                if nargin > 4
                    o.L = L;
                    
                    if nargin > 5
                        o.k = k;
                        
                        if nargin > 6
                            o.mass = mass;
                        
                            if nargin > 7
                                o.kin_h_function = kin_h_function;
                            end
                        end
                    end
                end
            end
            
            % == GUI ==
            o.x_plot = figure('Visible', 'OFF',...
                              'CloseRequestFcn', @(~,~) delete(o),...
                              'Position', [100, 200, 560, 420]);
            o.k_plot = figure('Visible', 'OFF',...
                              'CloseRequestFcn', @(~,~) delete(o),...
                              'Position', [660, 200, 560, 420]);
            o.slider = uicontrol(o.x_plot, 'Style', 'Slider',...
                                           'Units', 'normalized',...
                                           'Position', [0,-0.95,1,1]);
            
            if verLessThan('matlab', '7.13.0')
                event_name = 'ActionEvent';
            else
                event_name = 'ContinuousValueChange';
            end
            % Different Matlab versions use different event names.
            
            warning('off','MATLAB:addlistener:invalidEventName');
            addlistener(o.slider, event_name,@(~,~) o.internal_replot());
        end

        function set.V(o, value)
            o.V = value;
            o.recompute();
        end
        function set.L(o, value)
            o.L = value;
            o.recompute();
        end
        function set.k(o, value)
            o.k = value;
            o.redistribute();
        end
        function set.mass(o, value)
            o.mass = value;
            o.recompute();
        end
        
        function set.kin_h_function(o, value)
            o.kin_h_function = value;
            o.recompute();
        end
        
        function set.d(o, value)
            o.d = value;

            if isempty(o.d_mu) && isempty(o.d_sigma) %#ok<MCSUP>
                o.redistribute();
            end
        end
        function set.d_mu(o, value)
            o.d_mu = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_sigma) %#ok<MCSUP>
                o.redistribute();
            end
        end
        function set.d_sigma(o, value)
            o.d_sigma = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_sigma) %#ok<MCSUP>
                o.redistribute();
            end
        end
        
        function set.plot(o, value)
            o.plot = value;
            
            if value
                if(o.recompute_pending) %#ok<MCSUP>
                    o.internal_recompute();
                    o.internal_redistribute();
                elseif(o.redistribute_pending) %#ok<MCSUP>
                    o.internal_redistribute();
                end
                
                o.recompute_pending    = false; %#ok<MCSUP>
                o.redistribute_pending = false; %#ok<MCSUP>
                
                o.internal_replot();
                
                set(o.k_plot, 'Visible', 'ON'); %#ok<MCSUP>
                set(o.x_plot, 'Visible', 'ON'); %#ok<MCSUP>
            else
                set(o.k_plot, 'Visible', 'OFF'); %#ok<MCSUP>
                set(o.x_plot, 'Visible', 'OFF'); %#ok<MCSUP>
            end
        end
        function set.t_max(o, value)
            o.t_max = value;
            
            t_max_was_reset(o);
            % Since the response is different in both subclasses we have to
            % define there what to do!
        end
        
        function out = get.n(o)
            out = size(o.V,1);
        end

        function internal_recompute(~)
            error('Has to be implemented in a subclass!');
        end
        function internal_redistribute(~)
            error('Has to be implemented in a subclass!');
        end
        function internal_replot(~)
            error('Has to be implemented in a subclass!');
        end
        function t_max_was_reset(~,~)
            error('Has to be implemented in a subclass!');
        end
                
        % Create a new EVSimulation from this.
        function new = new_ev(o)
            new = EVSimulation(o.V, o.t_max, o.d_mu, o.d_sigma, o.L, o.k,...
                               o.mass, o.kin_h_function);
            new.d = o.d;
        end
        
        % Create a new StepwiseSimulation from this.
        function new = new_stepwise(o)
            new = StepwiseSimulation(o.V, o.t_max, o.d_mu, o.d_sigma, o.L, o.k,...
                                     o.mass, o.kin_h_function);
            new.d = o.d;
        end
        
        function delete(o)
            delete(o.x_plot);
            delete(o.k_plot);
        end
    end
    
    methods (Access = protected)
        function recompute(o)
            if o.plot
                o.internal_recompute();
                o.internal_redistribute();
                o.internal_replot();
            else
                o.recompute_pending = true;
            end
        end
        function redistribute(o)
            if o.plot
                o.internal_redistribute();
                o.internal_replot();
            else
                o.redistribute_pending = true;
            end
        end
        
        function replot(o)
            if o.plot
                o.internal_replot();
            end
        end
        
        function internal_recompute_kinetic_hamiltonian(o)
            % The kinetic hamiltonian:
            o.T = o.kin_h_function(o.n, o.L, o.mass);
        end
        function hamiltonian = get_hamiltonian(o, v_vec)
            assert(size(v_vec,2) == 1, 'Potential v_vec must be a column vector!');

            hamiltonian = diag(v_vec) + o.T;
        end
        function distribution = internal_redistribute_momentum(o)
            if ~isempty(o.d_mu) && ~isempty(o.d_sigma)
                o.d = sqrt(normpdf((1:o.n)', o.d_mu*o.n, o.d_sigma*o.n));
            end
            % Make sure d is set.
            
            if o.k ~= 0
                distribution = exp(1i*o.k*o.L*((1:o.n)'/o.n)) .* o.d(1:o.n);
            else
                distribution = o.d;
            end
            
            o.V_factor = max(abs(o.d))^2 / max(max(o.V));
        end
    end
    
    methods (Access = protected, Static = true)
        % Plot column vectors from 0 to 1:
        function plot_vec(varargin)
            n   = size(varargin{1},1);
            % Discretisation in x direction.
            len = length(varargin);
            % Length of arguments.
            
            x_values = (0:n-1)' / (n-1);
            % Go from 0 to 1 in n steps.

            plot_args = cell(2*len, 1);
            % Reserve memory for plot arguments.

            for i = 1:len
                plot_args{2*i - 1} = x_values;
                plot_args{2*i}     = varargin{i};
            end
            
            plot(plot_args{:}); %#ok<CPROP,PROP>
        end
        
        function plot_distr(distr, varargin)
            Simulation.plot_vec(abs(distr).^2, real(distr).^2, imag(distr).^2,...
                ... % Plot the sqares of the absolute, real and imaginary parts.
                                varargin{:});
                % ...And the rest of the arguments. (Like a potential.)
        end
    end
end

