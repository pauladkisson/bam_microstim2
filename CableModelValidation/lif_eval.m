function is_AP = lif_eval(Vmir)
    EL = -70e-3; %mV
    Vs = -50e-3; %mV
    if Vmir + EL >= Vs - 3e-3 %correction for different thresholds
        is_AP = 1;
    else
        is_AP = 0;
    end
end