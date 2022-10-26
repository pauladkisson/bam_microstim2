%%% Paul Adkisson
%%% 10/26/22
%%% Purpose: Validation of Extracellular Current Stimulation on LIF neurons
%%% using the Cable Equation model as a benchmark

%% Define Parameters
rattay_z_constants(1000*10^(-4)) %overwrite z value
load("rattay_constants.mat");
tstep = 0.0025; %ms
tmax = 5; %ms
I_el = -10; %uA
phase_dur = 0.3; %ms
tstart_ext = 1; %ms
tflip_ext = tstart_ext + phase_dur; %ms
tend_ext = tflip_ext + phase_dur; %ms
tstart_int = 0; %ms
tend_int = tmax;
V0 = V_rest;
y_ss = [V0, n_inf(V0), m_inf(V0), h_inf(V0)];
y0 = zeros(N, 4);
y0 = y0 + y_ss;
I_int = 0 * ones(N, 1) * 1e-6; %uA --> pA

%% Run Initial Test
tic;
[t, y, V_e] = rattayrun(y0, I_el, I_int, tstart_int, tend_int, ...
                           tstart_ext, tflip_ext, tend_ext, tmax, tstep, true);
toc
%% Plot Vm%%
nplot = 0:1:10;
xidx = 1:length(n);
xidx = xidx(any(n == nplot, 2));
cmap = turbo(length(nplot));
lin_cmap = turbo(1000);
fig = figure();
set(gca, 'Fontsize', 20)
colormap(lin_cmap);
hold on
for i = 1:length(nplot)
    idx = xidx(i);
    yplot = reshape(y(:, idx, 1), length(t), 1) - V_rest;
    plot(t, yplot, 'Color', cmap(i, :), 'Linewidth', 4)
end
hold off
xlabel("Time (ms)")
ylabel(["Membrane Polarization", "Relative to V_{rest} (mV)"])
title(["-10uA External Current Step Response", "(1mm away)"])
colorbar('Direction', 'reverse');

%% Hard-coded Voltage Step in LIF
C = C_m;
gL = g_l*pi*d*L;
EL = -70;
dt = tstep;
t = 0:dt:tmax;
trat = 0:tstep:tmax;
Vmir = (mean(V_e) - V_e(n==0));
I_inj = Vmir*gL;
V_step = ones(length(t), 1)*EL;
for i = 1:length(t)-1
    if t(i)<tstart_ext
        I_inj = 0;
    else
        if t(i) < tflip_ext
            I_inj = Vmir*gL;
            if t(i) == tstart_ext
                V_step(i) = V_step(i) + Vmir;
            end
        elseif t(i) <= tend_ext - tstep
            I_inj = -Vmir*gL;
            if t(i) == tflip_ext
                V_step(i) = V_step(i) - 2*Vmir;
            end
        else
            I_inj = 0;
            if abs(t(i)-tend_ext) < tstep/2
                V_step(i) = V_step(i) + Vmir;
            end
        end
    end
    if t(i)<tstart_int
        lif_I_int = 0;
    else
        lif_I_int = I_int(n==0);
    end
    dvdt = 1/C*( -gL*(V_step(i)-EL) + I_inj + lif_I_int);
    V_step(i+1) = V_step(i) + dvdt*dt;
end

%% Plot LIF step vs Cable Eq External Current Biphasic Pulse Response
Vrattay = y(:, n==0, 1) - V_rest;
figure;
set(gca, 'Fontsize', 20)
hold on 
plot(trat, Vrattay, "k-", 'LineWidth', 4)
plot(t, (V_step-EL), 'r--', 'LineWidth', 4)
xlabel("Time (ms)")
ylabel(["Membrane Polarization", "Relative to V_{rest} (mV)"])
legend(["Cable Equation", "LIF-step"])

lif_error = Vrattay - (V_step-EL);
lif_perc_error = lif_error / Vrattay;
figure;
set(gca, 'Fontsize', 20)
hold on
plot(t, lif_error, 'k-', 'Linewidth', 4)
xlabel("Time (ms)")
ylabel("Error (mV)")

%% CableEq: Depolarization vs Threshold at 5uA, 10uA, and 20uA
load("rattay_constants.mat");
show_prog = false;
z_bounds = [10, 20].*10^(-4);
z_thresh = z_bounds(1);
z_res = 0.1*10^(-4);
I_els = [-5, -10, -20]; %uA
I_ints = [0, 0.5, 1, 1.5, 2, 2.2, 2.4]*1e-6; %uA-->pA
depols = zeros(size(I_ints));
depol_thresholds = zeros(length(I_els), length(I_ints));
for j = 1:length(I_els)
    I_el = I_els(j);
    fprintf("External Current=%0.0fuA \n", I_el);
    for i = 1:length(I_ints)
        I_int = I_ints(i);
        fprintf("Internal Current=%0.1fpA \n", I_int*1e6)
        y = z_rattay(100*1e-4, 0, I_int, true);
        depol = max(y(:, 1, 1)) - V_rest;
        depols(i) = depol;
        fprintf("Depolarization=%0.0fmV \n", depol)
        z_rat = @(z) z_rattay(z, I_el, I_int, show_prog);
        try
            %new threshold is greater than previous to avoid numerical instability
            z_bounds(1) = depol_thresholds(j, i-1);
        catch %first I_int for a new I_el
            assert(i==1, "Some other error has occured")
            try
                z_bounds(1) = depol_thresholds(j-1, 1);
            catch %first I_el
                assert(j==1, "Some other error has occured")
            end
        end
        z_bounds(2) = z_bounds(1)*2; %new threshold is no more than double the previous one
        z_thresh = bisect_search(z_rat, @rattay_eval, z_bounds, z_res);
        fprintf("Threshold Distance = %0.1fum \n", z_thresh*1e4)
        depol_thresholds(j, i) = z_thresh;
    end
end
save("threshold_correction.mat", 'depol_thresholds', 'depols', 'I_ints', 'I_els')

%% Optimize threshold correction factor
load("threshold_correction.mat") 
bounds = [0, 1]; %thresh correction bounds
options = optimset('Display','iter');
eval_fn = @(thresh_cor) eval_thresh_cor(I_els, depol_thresholds, thresh_cor);
thresh_cor_star = fminbnd(eval_fn, bounds(1), bounds(2), options);
save("threshold_correction.mat", 'depol_thresholds', 'depols', 'I_ints', ...
    'I_els', 'thresh_cor_star');

%% LIF Approximation: Depolarization vs Threshold at 5uA, 10uA, and 20uA
load("threshold_correction.mat");
C = 0.5*1e-9; %nF
gL = 25*1e-9; %nS
EL = -70e-3; %mV
z_bounds = [0, 400]*1e-4;
z_res = 1e-4;
lif_depol_thresholds = zeros(length(I_els), length(depols));
for j = 1:length(I_els)
    I_el = I_els(j);
    fprintf("External Current Amplitude = %0.0fuA \n", I_el)
    for i = 1:length(depols)
        depol = depols(i);
        fprintf("Depolarization=%0.0fmV \n", depol)
        mirror_est = @(z) mir_est(z, I_el, thresh_cor_star, depol);
        z_thresh = bisect_search(mirror_est, @lif_eval, z_bounds, z_res);
        fprintf("Threshold Distance = %0.1fum \n", z_thresh*1e4)
        lif_depol_thresholds(j, i) = z_thresh;
    end
end

%% Plot Goodness of Fit
figure;
set(gca, 'Fontsize', 20)
hold on
lin_I_els = 0:-0.01:I_els(end);
cmap = turbo(length(lin_I_els));
for j = 1:length(I_els)
    I_el = I_els(j);
    I_color = cmap(lin_I_els==I_el, :);
    plot(depols(1:end-1), depol_thresholds(j, :)*10^4, 'o-', 'Linewidth', 4, 'Color', I_color)
    plot(depols(1:end-1), lif_depol_thresholds(j, :)*10^4, '--', 'Linewidth', 4, 'Color', I_color)
end
colormap(cmap)
c = colorbar('Ticks', [0, I_els./I_els(end)], 'TickLabels', [0, compose("%0.0f", abs(I_els))]);
c.Label.String = 'Pulse Amplitude (uA)';
xlabel("Depolarization (mV)")
ylabel("Activation Threshold (um)")
legend(["Cable Eq", "LIF"])
title(sprintf("Optimal Threshold Correction = %0.03f", thresh_cor_star))

figure;
set(gca, 'Fontsize', 20)
hold on
for j = 1:length(I_els)
    I_el = I_els(j);
    I_color = cmap(lin_I_els==I_el, :);
    err = depol_thresholds(j, :) - lif_depol_thresholds(j, :);
    plot(depols(1:end-1), err*1e4, 'o-', 'Linewidth', 4, 'Color', I_color)
end
colormap(cmap)
c = colorbar('Ticks', [0, I_els./I_els(end)], 'TickLabels', [0, compose("%0.0f", abs(I_els))]);
c.Label.String = 'Pulse Amplitude (uA)';
xlabel("Depolarization (mV)")
ylabel("Absolute Error (um)")
title(sprintf("Optimal Threshold Correction = %0.03f", thresh_cor_star))

figure;
set(gca, 'Fontsize', 20)
hold on
for j = 1:length(I_els)
    I_el = I_els(j);
    I_color = cmap(lin_I_els==I_el, :);
    err = depol_thresholds(j, :) - lif_depol_thresholds(j, :);
    rel_err = err ./ depol_thresholds(j, :);
    plot(depols, rel_err*100, 'o-', 'Linewidth', 4, 'Color', I_color)
end
colormap(cmap)
c = colorbar('Ticks', [0, I_els./I_els(end)], 'TickLabels', [0, compose("%0.0f", abs(I_els))]);
c.Label.String = 'Pulse Amplitude (uA)';
xlabel("Depolarization (mV)")
ylabel("Relative Error (%)")
title(sprintf("Optimal Threshold Correction = %0.03f", thresh_cor_star))
