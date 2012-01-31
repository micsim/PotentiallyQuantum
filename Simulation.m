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
        V_factor;
        
        x_plot;
        k_plot;
        slider;
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
                              'CloseRequestFcn', @(~,~) delete(o));
            o.k_plot = figure('Visible', 'OFF',...
                              'CloseRequestFcn', @(~,~) delete(o));
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

            if isempty(o.d_mu) && isempty(o.d_sigma)
                o.redistribute();
            end
        end
        function set.d_mu(o, value)
            o.d_mu = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_sigma)
                o.redistribute();
            end
        end
        function set.d_sigma(o, value)
            o.d_sigma = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_sigma)
                o.redistribute();
            end
        end
        
        function set.plot(o, value)
            o.plot = value;
            
            if value
                o.internal_replot();
                
                set(o.x_plot, 'Visible', 'ON');
                set(o.k_plot, 'Visible', 'ON');
            else
                set(o.x_plot, 'Visible', 'OFF');
                set(o.k_plot, 'Visible', 'OFF');
            end
        end
        function set.t_max(o, value)
            o.t_max = value;
            
            o.replot();
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
        
        function delete(o)
            delete(o.x_plot);
            delete(o.k_plot);
        end

    end
    
    methods (Access = protected)
        % Redo all computations.
        function set_property_recompute(o, name, value)
            o.(name) = value;
            
            o.recompute();
        end
        
        % Redo the evolution including changes from the starting k
        function set_property_reevolve(o, name, value)
            o.(name) = value;
            
            o.reevolve();
        end
        
        function recompute(o)
            if o.plot
                o.internal_recompute();
                o.internal_redistribute();
                o.internal_replot();
            end
        end
        function redistribute(o)
            if o.plot
                o.internal_redistribute();
                o.internal_replot();
            end
        end
        
        function replot(o)
            if o.plot
                o.internal_replot();
            end
        end
    end
end
