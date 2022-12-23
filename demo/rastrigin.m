function [y] = rastrigin(x)
    if isrow(x)
        d = size(x,2);
    else
        d = size(x,1);
    end
    sq = x.^2;   
    y = 10*d + sum(sq - 10*cos(2*pi*x));
end