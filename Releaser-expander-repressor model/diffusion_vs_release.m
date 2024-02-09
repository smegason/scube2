% generates comparison of diffusion only and diffusion+release models
% as the production rate is varied

phis = 10 .^ (-1 : 0.02 : 3);

lengths = [200, 300];

load('variables\self_enhanced.mat');
handle = @models.HS_self_enhanced;

% set up expansion only and both parameter sets
paramsModel_both = paramsModel;
paramsModel.Phi_HS_rel = 1;
paramsModel_expansion = paramsModel;

paramsSys.t_min = 1440;
paramsSys.max_t = 1440;

data = zeros(2, numel(phis));

for i = 1 : numel(phis)
    data(1, i) = get_scaling(paramsModel_expansion, paramsSys, handle, ...
        [1, 1] * phis(i), lengths);
    data(2, i) = get_scaling(paramsModel_both, paramsSys, handle, ...
        [1, 1] * phis(i), lengths); 
end

plot(phis, data.');
ax = gca;
ax.XScale = 'log';
xlabel('HH production rate');   
ylabel('Scaling');
legend('Diffusion only', 'Diffusion and Release');

function eta = get_scaling(paramsModel, paramsSys, handle, phis, lengths)
    % runs model for the two given lengths and computes the scaling metric
    % between the two steady state gradients
    
    paramsSys.L = lengths(1);
    paramsModel.a_H = phis(1);
    model = handle(paramsModel, paramsSys);
    [~, y_short] = run_model(paramsSys, model);
    A_short = squeeze(y_short(:, 1, :));
    
    paramsSys.L = lengths(2);
    paramsModel.a_H = phis(2);
    model = handle(paramsModel, paramsSys);
    [~, y_long] = run_model(paramsSys, model);
    A_long = squeeze(y_long(:, 1, :));

    eta = metrics.L2_area(A_short, A_long);
end
        
