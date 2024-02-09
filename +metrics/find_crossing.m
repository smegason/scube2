function n = find_crossing(y, threshold)
    % given state vector y, finds the first cell n such that y(n) is less
    % than the threshold; assumes y is monotically decreasing
    % interpolates between cells; returns an n that is not neccesarily a
    % whole number
    
    n_max = find(y < threshold, 1, 'first');
    
    if isempty(n_max)
        n = Inf; % never crossed threshold
        return
    elseif n_max == 1
        n = 0;  % crossed threshold before beginning of y
        return
    end
        
    % get changes across cell n and n - 1
    Dy = y(n_max - 1) - y(n_max);
    dy = threshold - y(n_max);
    
    n = n_max - (dy/Dy) - 1/2;
    
end