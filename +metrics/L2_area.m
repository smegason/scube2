% computes \int_0^L (f - g)^2 dx between the two curves, then 
% normalizes to \int_0^L f^2 (x) dx

function eta = L2_area(y1, y2)
    difference = sum((y1(end, :) - y2(end, :)) .^ 2);
    normalization = sum(y1(end, :) .^2);
    eta = difference / normalization;
end
