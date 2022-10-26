function is_AP = rattay_eval(y)
    if max(y(:, end, 1)) > 0 || any(isnan(y(:, end, 1)))
        is_AP = 1;
    else
        is_AP = 0;
    end
end