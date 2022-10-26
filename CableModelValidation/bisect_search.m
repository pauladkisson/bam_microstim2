%%% Implement Bisection Search
%%% Inputs:
%%%     my_func : function to be evaluated
%%%     eval : evaluates the output of @my_func and returns True if x achieves
%%%         desirable outcome, otherwise False
%%%     bounds : [lower bound, upperbound]
%%%     resolution : change in x defining convergence
%%% Output:
%%%     x : threshold
function x = bisect_search(my_func, eval, bounds, resolution)
    x_low = bounds(1);
    x_high = bounds(2);
    x_test = x_low + (x_high - x_low)/2;
    while abs(x_test-x_low) > resolution
        if eval(my_func(x_test))
            x_low = x_test;
        else
            x_high = x_test;
        end
        x_test = x_low + (x_high - x_low)/2;
    end
    x = x_test;
end