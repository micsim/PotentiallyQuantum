function plot_timedep_slider(psi, n_t, dt, v_mat)
    fig = figure();
    slider = uicontrol(fig, 'Style', 'Slider', 'Units', 'normalized',...
        'Position', [0,-0.95,1,1]);
    warning('off','MATLAB:addlistener:invalidEventName');
    addlistener(slider,'ActionEvent',@(~,~) slider_moved());
    slider_moved();

    function slider_moved()
        value = floor(get(slider, 'Value')*(n_t-1) + 1);

        plot_timedep(psi(:,value), v_mat(:,value));
    end
end
