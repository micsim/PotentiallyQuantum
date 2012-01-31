function plot_slider(e, E, d, v, t_max, K)
    x_plot = figure();
    slider = uicontrol(x_plot, 'Style', 'Slider', 'Units', 'normalized',...
        'Position', [0,-0.95,1,1]);
    warning('off','MATLAB:addlistener:invalidEventName');
    
    if verLessThan('matlab', '7.13.0')
        event_name = 'ActionEvent';
    else
        event_name = 'ContinuousValueChange';
    end
    % Different Matlab versions use different event names.
    addlistener(slider, event_name,@(~,~) slider_moved());

    plot_distr(E, d, v);

    has_K = (nargin == 6);

    if has_K
        k_plot = figure();
        plot_distr(K, d, 0*v);

        k_plot2 = figure();
        plot_timedep(fft(E*d), 0*v);
    end

    function slider_moved()
        value = get(slider, 'Value')*t_max;
        E_new = evolve_eigenvectors(e, E, value);

        set(0,'CurrentFigure',x_plot);
        plot_distr(E_new, d, v);

        if has_K
            K_new = evolve_eigenvectors(e, K, value);

            set(0,'CurrentFigure',k_plot);
            plot_distr(K_new, d, 0*v);
            set(0,'CurrentFigure',k_plot2);
            plot_timedep(fft(E_new*d), 0*v);
        end
    end
end
