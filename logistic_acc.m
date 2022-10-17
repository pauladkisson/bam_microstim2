%%% Implement Logistic Function
%%% Paul Adkisson
%%% 1/17/22

function w = logistic_acc(coeffs, c)
    a = coeffs(1);
    b = coeffs(2);
    X = a + b*c;
    w = 1 ./ (1 + exp(-X));
end