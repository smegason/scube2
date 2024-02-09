% contains code that, given a model class, runs the model

function [t, y, out_flag] = run_model(paramsSys, model, preequilibrium, set_concentration)
    % inputs: paramsSys (struct), model (model object), preeq (bool)
    % given a particular model class, loads and runs the model
    % if preeq is true, runs a preequilibrium step; otherwise
    % just uses the 0D initial condition as the initial condition
    % if set_concentration is true, sets the hedgehog concentration
    % at the origin to model.H_0 at the beginning of the simulation
    
    n = model.numSpecies;
    N = model.N;
    
    if nargin == 2
        preequilibrium = false;  % default false
        set_concentration = false;
    end
    
    if preequilibrium
        pre_ode_options = odeset('NonNegative', 1 : n-1);
        pre_ICs = model.get_initial_conditions_0D();
        [~, pre_y] = ode15sTillSs(paramsSys, @model.model_0D, pre_ICs, ...
            model.get_Ss_species_0D(), pre_ode_options);
        ICs = ones(N, 1) * pre_y;
        ICs = reshape(ICs, [], 1);
    else
        ICs = ones(N, 1) * model.get_initial_conditions_0D();
        ICs = reshape(ICs, [], 1);
    end
    
    if set_concentration
        ICs(1) = model.H_0;
    end

    % set up ODE-solvers
    ode_options = odeset('JPattern', model.get_jacobian_sparsity(), 'NonNegative', 1 : N * n);
    
    % run 1D model
    [t, y, out_flag] = ode15sTillSs(paramsSys, @model.model_1D, ICs, ...
        model.get_Ss_species_1D(), ode_options);
    
    % reshape y into "matrix" form
    y = reshape(y, length(t), model.numSpecies, model.N);
end

function [t, y, outFlag] = ode15sTillSs(paramsSys, modelHandle, y0, species2Equil, options)
% uses ode15s solver to solve ODE until it reaches steady state
% divides tspan into chunks and checks every check_interval to see if 
% steady state is reached
% outFlag
%   1 -- reached steady state
%   2 -- did not reach steady state and hit maxTime
%   3 -- error was generated    

    t_span = 0 : paramsSys.t_step : paramsSys.max_t;
    check_int = paramsSys.sS_check_interval;

    % set up variables
    y0Temp = y0; % initial conditions
    t = 0;
    y(1, :) = y0;
    outFlag = 2;
   
    for i = 1 : ceil((length(t_span) - 1) / check_int)
        % Run model for time t_chunk
        try
            t_chunk = t_span((i - 1) * check_int + 1 : min(length(t_span), i * check_int + 1));
            [tTemp, yTemp] = ode15s(modelHandle, t_chunk, y0Temp, options);
            t = cat(1, t, tTemp(2:end)); % add tTemp to end of t{i}
            y = cat(1, y, yTemp(2:end,:)); % add yTemp to the end of y{i}
        catch exception
            exception.getReport
            outFlag = 3;
        end
        
        % Check for a break condition       
        if outFlag == 3
            break
        end
        
        if isOutputSteadyState(t, y, species2Equil, paramsSys)
            outFlag = 1;
            break
        end
        
        % otherwise proceed to next time chunk
        y0Temp = y(end,:);
    end
end

function flag = isOutputSteadyState(t, y, species2Equil, paramsSys)
% if (1) we've exceeded t_min and the number of lookbacksteps and
% (2) the *maximum* fractional rate change of every species2Equail is
% below sSThreshold, returns true; false otherwise

    lookback_steps = paramsSys.sS_lookbackSteps;
    num_steps = size(y, 1);

    if t(end) < paramsSys.t_min || num_steps < lookback_steps + 1
        flag = 0;
        return
    end
    
    lookback = y(num_steps - lookback_steps : num_steps, :);
    frac_diffs = lookback(:, species2Equil) ./ lookback(1, species2Equil) - 1;
    frac_derivs = frac_diffs / paramsSys.t_step;
    
    if max(frac_derivs, [], 'all') < paramsSys.sS_threshold
        flag = 1;
    else
        flag = 0;
    end
end
