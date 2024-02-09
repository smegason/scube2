% runs a model for two different lengths of the patterning
% field, the plots the steady-state perception against either
% position or relative position

t_end = 1440;


%% Uniform degradation

lengths = [200, 300];
phis = [1, 1];
handle = @models.HS;

load('variables/uniform_degradation');
make_plots(paramsModel, paramsSys, handle, phis, lengths, t_end, 1)

% now, do release only
paramsModel.D_HS_rel = 1;
make_plots(paramsModel, paramsSys, handle, phis, lengths, t_end, 3)

% now, do zero scaling
paramsModel.Phi_HS_rel = 1;
make_plots(paramsModel, paramsSys, handle, phis, lengths, t_end, 5)

%%
function make_plots(paramsModel, paramsSys, handle, phis, lengths, t_end, fig_num)

    paramsSys.t_min = t_end;
    paramsSys.max_t = t_end;

    disp('Running model for short field...')
    [A_short, t_short] = get_dynamics(1, paramsModel, paramsSys, handle, phis(1), lengths(1));
    disp('Running model for long field...')
    [A_long, t_long] = get_dynamics(1, paramsModel, paramsSys, handle, phis(2), lengths(2));

    % as time, pick the longest 
    if length(t_long) > length(t_short)
        t = t_long;
    else
        t = t_short;
    end

    figure(fig_num)
    plot_data(paramsSys, lengths, A_short, A_long, 0, t, length(t));
    figure(fig_num + 1)
    plot_data(paramsSys, lengths, A_short, A_long, 1, t, length(t));
    
    scaling = metrics.domains(A_long, A_short);
    disp(['Scaling score: ' num2str(scaling)]);

end


%%

function [output, t] = get_dynamics(index, paramsModel, paramsSys, handle, phi, L)
    paramsSys.L = L;
    paramsModel.a_H = phi;
    model = handle(paramsModel, paramsSys);
    [t, y] = run_model(paramsSys, model);
    
    A = zeros(length(t), paramsSys.N);
    for i = 1 : length(t)
        A(i, :) = model.get_perception(squeeze(y(i, :, :)));
    end

    if strcmp(index(1), 'p')
        output = A;
    else
        output = squeeze(y(:, index, :));
    end
end


function plot_data(paramsSys, lengths, A_short, A_long, scale_axis, t, i);
    % plots the ith timestepfor all data sets in the cell array 
    % data at the ith timstep
   
    % set up axes
    y_min = min(min(A_short, [], 'all'), min(A_long, [], 'all'));
    y_max = max(max(A_short, [], 'all'), max(A_long, [], 'all'));
    a = 0.05;  % how much extra space to show
    
    y_min = max(y_min, -7);  % so it doesn't break when on a log scale
    % y_max = 2;
    
    % y_axis = [y_min - a*(y_max - y_min), y_max + a*(y_max - y_min)];
    y_axis = [y_min, y_max + a*(y_max - y_min)];
    x_axis = lengths(2) / paramsSys.N * (1 : paramsSys.N);
    
    clf;
    hold on;
    
    gradient_short = squeeze(A_short(i, :));
    gradient_long = squeeze(A_long(i, :));
    
    if scale_axis
        plot(x_axis / lengths(2), gradient_short);
        plot(x_axis / lengths(2), gradient_long);
        xlabel('Relative position x/L');
    else
        x_axis_short = 1 : ceil(lengths(1) / lengths(2)) * paramsSys.N;
        interped = interp1((1 : paramsSys.N) * (lengths(1) / paramsSys.N), ...
            gradient_short, x_axis_short);
        plot(x_axis_short, interped);
        
        plot(x_axis, gradient_long);
        xlabel('Position (um)');
    end
    
    ylim(y_axis);
    ylabel('Scaled ligand concentration');
    hold off;
end
