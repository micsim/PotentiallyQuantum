function plot_timedep_slider(psi, n_t, dt, v_mat)
    x_plot = figure();
    k_plot = figure();
    slider = uicontrol(x_plot, 'Style', 'Slider', 'Units', 'normalized',...
        'Position', [0,-0.95,1,1]);
    warning('off','MATLAB:addlistener:invalidEventName');
    addlistener(slider,'ActionEvent',@(~,~) slider_moved());
    slider_moved();

    function slider_moved()
        value = floor(get(slider, 'Value')*(n_t-1) + 1);

        set(0,'CurrentFigure',x_plot);
        plot_timedep(psi(:,value), v_mat(:,value));
        set(0,'CurrentFigure',k_plot);
        plot_timedep(fft(psi(:,value)), 0*v_mat(:,1));
    end
end
