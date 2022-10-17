%%% Implement Weibull Function
%%% Paul Adkisson
%%% 9/15/21

function w = weibull(coeffs, c)
    a = coeffs(1);
    b = coeffs(2);
    w = 1 - 0.5*exp(-(c/a).^b);
end