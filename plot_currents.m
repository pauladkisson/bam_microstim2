%% Simulation Parameters
clear;
sim_name = "record_currents";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))

% Adjust stim_amps
pulse_amps = [-10*1e-6];
dc_amps = [-1.4, 0, 1.4]*1e-6;
stim_amps = [pulse_amps, dc_amps];default_colors = get(gca, "colororder");
stim_coherences = [-55, -55, 0, 30] ./ 100;


%% Assess Different stim conditions side by side
I_nets = zeros(length(t), length(stim_amps));
stim_frs = zeros(length(t), length(stim_amps));
for j = 1:length(stim_amps)
    fprintf("j = %0.0f \n", j)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
    if pulse
        trialpath = strcat(sim_path, sprintf("/data/%0.2fuA_pulse/c=%0.3f/trial2.mat", ...
            [stim_amp*1e6, stim_coherences(j)]));
    else
        trialpath = strcat(sim_path, sprintf("/data/%0.2fuA_galvanic/c=%0.3f/trial2.mat", ...
            [stim_amp*1e6, stim_coherences(j)]));
    end
    load(trialpath)
    % Filter the currents in time
    avg_win_size = 50e-3;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    I_AMPA_ext = filter(w, 1, I_AMPA_ext);
    I_AMPA_rec = filter(w, 1, I_AMPA_rec);
    I_NMDA = filter(w, 1, I_NMDA);
    I_GABA = filter(w, 1, abs(I_GABA));
    I_LEAK = filter(w, 1, abs(I_LEAK));
    I_rec_exc = I_AMPA_rec + I_NMDA;
    I_net = I_rec_exc - I_GABA;
    I_nets(:, j) = I_net;
    
    keys = 1:size(Vm, 2);
    t_idx = 1:size(Vm, 1);
    recspikes = containers.Map;
    for key = keys
        recspikes(int2str(key)) = t_idx(Vm(:, key)==0);
    end
    [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
    stim_frs(:, j) = pop_frs(:, 1);
end

%% plot
task_mask = t >= 0 & t < 5;

figure;
hold on
plot(t(task_mask), stim_frs(task_mask, 1), 'Color', default_colors(7, :), 'Linewidth', 2)
plot(t(task_mask), stim_frs(task_mask, 2), 'Color', default_colors(5, :), 'Linewidth', 2)
plot(t(task_mask), stim_frs(task_mask, 3), 'Color', 'k', 'Linewidth', 2)
plot(t(task_mask), stim_frs(task_mask, 4),  'Color', default_colors(6, :), 'Linewidth', 2)
xlabel("Time (s)")
ylabel("P1 Firing Rate (spk/s)")
title("P1 Loses (Trial 2)")
set(gca, 'Fontsize', 20)

figure;
hold on
plot(t(task_mask), I_nets(task_mask, 1)*1e12, 'Color', default_colors(7, :), 'Linewidth', 2)
plot(t(task_mask), I_nets(task_mask, 2)*1e12, 'Color', default_colors(5, :), 'Linewidth', 2)
plot(t(task_mask), I_nets(task_mask, 3)*1e12, 'Color', 'k', 'Linewidth', 2)
plot(t(task_mask), I_nets(task_mask, 4)*1e12,  'Color', default_colors(6, :), 'Linewidth', 2)
yline(0, 'k--')
xlabel("Time (s)")
ylabel("Network Current (pA)")
title("P1 Loses (Trial 2)")
set(gca, 'Fontsize', 20)


