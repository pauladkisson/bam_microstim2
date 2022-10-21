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
dc_amps = [-1, 0]*1e-6;
stim_amps = [pulse_amps, dc_amps];
%}

%% Plot Firing Rates
%ex_c = -25.6/100;
ex_c = 0;
ex_trial = 1;
ex_stim_j = 3;
plot_name = "single_stim"; % 'single_stim' or 'subplot' or 'p1_only'
plot_frs(sim_name, pulse_amps, stim_amps, p, t, t_task, t_taskoff, default_colors, ...
                  ex_stim_j, ex_c, ex_trial, plot_name);


%%
top_N = num_group*0.1;
ex_neurons = [5, 6];
ex_stim_j = 1;
ex_brain = 1;
ex_trial = 1;
plot_name = "grouped_stim"; % or 'single_stim' or 'grouped_stim'
plot_rasters(sim_name, pulse_amps, stim_amps, ex_neurons, t, t_task, t_taskoff, stim_freq, default_colors, top_N, ex_stim_j, ex_brain, ex_c, ex_trial, plot_name)
%}

%{
win_size = floor(0.250 / dt); %250ms moving window
cv_window = t >= 2.5 & t<3; %Plotting window
%cv_window = t >= t_task & t<t_taskoff; %Plotting window
ex_neuron = 7;
ex_brain = 1;
ex_trial = 1;
top_N = floor(num_group);
plot_name = "p1_wins"; % or 'ex_trial' or 'p1_wins'
sim_names = ["FeedforwardGamma", "FeedforwardGamma"];
%sim_names = ["%Activation Equivalence Connected", "%Activation Equivalence Disconnected"];
plot_cv(sim_name, sim_names, pulse_amps, stim_amps, t, N, top_N, num_group, ...
                 win_size, cv_window, default_colors, ex_brain, ex_c, ex_trial, ...
                 ex_neuron, brains, num_brains, pulse_coherences, galvanic_coherences, control_coherences, ...
                 start_trial, end_trial, num_trials, plot_name)
%}

%{
sim_names = ["FeedForwardGamma", "FeedforwardGamma"];
idx_diff = stim_ind+1;% how far off timing is from pulse timing + 1 to account for t(1) = 0
plot_phaselock(sim_names, pulse_amps, stim_amps, t, t_task, t_taskoff, stim_freq, num_group, ...
                        idx_diff, default_colors, brains, num_brains, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials)
%}

%{
N_start = 1;
N_end = floor(num_group);
win_start = 2.5;
win_stop = 3;
c_win = 300*1e-6;
c = 0;
sim_names = ["thresh_cor=0.211 Con", "thresh_cor=0.211 Discon"];
plot_sync(sim_names, pulse_amps, stim_amps, t, num_group, ...
                        brains, num_brains, N_start, N_end, ...
                        win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials, default_colors)
%}


win_start = t_task + stim_ind*dt; % to account for onset spike of pulse
win_stop = t_task + 0.1;
ex_c = 0;
%  plot_name = 'ex_c' or 'p1_wins' or 'p1_loses'
plot_name = "ex_c";
num_brains = length(brains);
plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, num_affected, ...
                     win_start, win_stop, default_colors, num_brains, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     start_trial, end_trial, num_trials, plot_name);
%}

%{
win_start = 2.5;
win_stop = 3;
ex_c = 0;
f0 = 40;
%  plot_name = 'ex_c' or 'p1_wins' or 'p1_loses'
plot_name = "p1_wins";
num_brains = length(brains);
plot_oscillation(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, ...
                        win_start, win_stop, default_colors, brains, num_brains, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        f0, start_trial, end_trial, num_trials, plot_name);
%}

%{
plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, brains, ...
               num_brains, num_batch, ...
               pulse_coherences, galvanic_coherences, control_coherences)
%}
