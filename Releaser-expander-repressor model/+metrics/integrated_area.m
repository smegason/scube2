% computes \int_t_start^t_end \int_0^L (f - g)^2 dx dt between the two curves, then 
% normalizes to \int t_start^t_end \int_0^L f^2 (x) dx dt

function eta = integrated_area(y1, y2, t_start, t_step, t_end)

    if nargin == 2  % defaults
        t_start = 0;
        t_step = 1;
        t_end = max(height(y1), height(y2)) - 1;
    end
    
    start_ind = round(t_start / t_step) + 1;
    end_ind = round(t_end / t_step) + 1;
    
    % if the given y's aren't long enough, assume we've reached steady
    % state and just "extend" the y's downward
    if height(y1) < end_ind
        extension = ones(end_ind - height(y1), 1) * y1(end, :);
        y1 = [y1; extension];
    end
        
    if height(y2) < end_ind
        extension = ones(end_ind - height(y2), 1) * y2(end, :);
        y2 = [y2; extension];
    end
    
    y1 = y1(start_ind:end_ind, :);
    y2 = y2(start_ind:end_ind, :);
    
    difference = sum((y1 - y2) .^ 2, 'all');
    normalization = sum(y1 .^2, 'all');
    eta = difference / normalization;
end
