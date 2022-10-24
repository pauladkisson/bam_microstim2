%% Simulation Parameters
clear;
sim_name = "Brainless_m=0_Discon";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
figure('visible', 'off');
default_colors = get(gca, "colororder");

start_trial = 1;
end_trial = 28;
num_trials = end_trial - start_trial + 1;
num_batch = 3;

control_coherences = 0;
galvanic_coherences = 0;
pulse_coherences = 0;

pulse_amps = [-10*1e-6];
dc_amps = [-1.4, 0]*1e-6;
stim_amps = [pulse_amps, dc_amps];

%% Plot Firing Rates
ex_c = 0;
ex_trial = 1;
ex_stim_j = 2;
plot_name = "p1_only"; % 'single_stim' or 'subplot' or 'p1_only'
plot_frs(sim_name, pulse_amps, stim_amps, p, f, N, N_E, t, t_task,...
                  t_taskoff, default_colors, ex_stim_j, ex_c, ex_trial, plot_name);


%% Plot Rasters
top_N = num_group*0.1;
ex_neurons = [3];
tlim = [1, 3];
ex_stim_j = 1;
ex_brain = 1;
ex_trial = 1;
ex_c = 0;
plot_name = "grouped_stim"; %'ex_trial' or 'single_stim' or 'grouped_stim'
plot_rasters(sim_name, pulse_amps, stim_amps, ex_neurons, t, ...
                      t_task, t_taskoff, stim_freq, default_colors, top_N, ...
                      ex_stim_j, ex_c, ex_trial, tlim, plot_name);
                  
%% Plot Coefficeint of Variation (CV)
win_size = floor(0.250 / dt); %250ms moving window
cv_window = t >= 2.5 & t<3; %Plotting window
ex_neuron = 7;
ex_trial = 1;
top_N = floor(0.1*num_group);
plot_name = "p1_wins"; %ex_neuron or 'ex_trial' or 'p1_wins'
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
plot_cv(sim_name, sim_names, pulse_amps, stim_amps, t, N, top_N, num_group, ...
                 num_affected, win_size, cv_window, default_colors, ex_c, ex_trial, ...
                 ex_neuron, pulse_coherences, galvanic_coherences, control_coherences, ...
                 start_trial, end_trial, num_trials, plot_name);


%% Plot Phaselocking to Pulses
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
idx_diff = stim_ind+1;% how far off timing is from pulse timing + 1 to account for t(1) = 0
win_start = 2.5;
win_stop = 3;
top_N = floor(1*num_group);
plot_phaselock(sim_names, pulse_amps, stim_amps, t, t_task, t_taskoff, stim_freq, ...
                        num_group, num_affected, top_N, win_start, win_stop, idx_diff, ...
                        default_colors, start_trial, end_trial, num_trials, ...
                        pulse_coherences, galvanic_coherences, control_coherences);

%% Plot Synchrony
N_start = 1;
N_end = floor(num_group);
win_start = 2.5;
win_stop = 3;
c_win = 300*1e-6;
c = 0;
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
symmetric = true;
plot_sync(sim_names, pulse_amps, stim_amps, t, num_group, N_start, ...
                        N_end, win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials, symmetric);

%% Plot Firing Rate Distribution over distance from electrode
win_start = t_task + stim_ind*dt; % to account for onset spike of pulse
win_stop = t_task + 0.1;
ex_c = 0;
%  plot_name = 'ex_c' or 'p1_wins' or 'p1_loses'
plot_name = "ex_c";
plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     start_trial, end_trial, num_trials, plot_name);

%% Plot Decisions and Decision Times
plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, num_batch, ...
               pulse_coherences, galvanic_coherences, control_coherences, num_trials);
