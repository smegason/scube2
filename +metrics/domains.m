% makes n evenly sized domains from the first gradient, then
% computes the average shift in the position of each gradient

function eta = domains(y1, y2, n)
    if nargin == 2
        n = 4;  % by default
    end

    y1 = squeeze(y1(end, :));
    y2 = squeeze(y2(end, :));
    
    N = length(y1);
    sum_delta = 0;
    
    for i = 1 : n
        % generate threshold
        n_cell = round(N * i/(n + 1));
        T = y1(n_cell);
        n_shifted = find_crossing(y2, T);
        
        sum_delta = sum_delta + abs(n_cell - n_shifted) / N;
    end
    
    eta = sum_delta / n;
end


function n = find_crossing(y, threshold)
    % given state vector y, finds the first cell n such that y(n) is less
    % than the threshold; assumes y is monotically decreasing
    % interpolates between cells; returns an n that is not neccesarily a
    % whole number
    
    n_max = find(y < threshold, 1, 'first');
    
    if isempty(n_max)
        n = length(y) + 1; % never crossed threshold
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