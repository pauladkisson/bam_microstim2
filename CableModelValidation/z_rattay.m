function y = z_rattay(z, I_el, I_int, show_prog)
    rattay_z_constants(z) %overwrite z value
    load("rattay_constants.mat")
    fprintf("z=%0.1f \n", z*10^4)
    tstep = 0.0025; %ms
    tmax = 10; %ms
    phase_dur = 0.3; %ms
    tstart_ext = 5; %ms
    tflip_ext = tstart_ext + phase_dur; %ms
    tmax_ext = tflip_ext + phase_dur; %ms
    tstart_int = 0; %ms
    tmax_int = tmax;
    V0 = V_rest;
    y_ss = [V0, n_inf(V0), m_inf(V0), h_inf(V0)];
    y0 = zeros(N, 4);
    y0 = y0 + y_ss;
    [~, y, ~] = rattayrun(y0, I_el, I_int, tstart_int, tmax_int, ...
                           tstart_ext, tflip_ext, tmax_ext, tmax, tstep, show_prog);
end