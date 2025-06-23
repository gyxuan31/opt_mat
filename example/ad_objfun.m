function obj = objfun(x, m, n, P, R, B)
    D = reshape(x, [m, n]);
    s = R(:) .* sum(P .* D, 2);
    obj = -sum(s);
end

