%% Simulation Parameters
clear;
sim_name = "iScience_Con";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
figure('visible', 'off');
default_colors = get(gca, "colororder");

start_trial = 1;
end_trial = 56;
num_trials = end_trial - start_trial + 1;
num_batch = 7;

%control_coherences = 0;
%galvanic_coherences = 0;
%pulse_coherences = 0;
%anodic_coherences = 0;

%{
control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, ...
                      3.2, 6.4, 12.8, 25.6, 51.2] / 100;
pulse_coherences = [-100, -65, -60, -59, -58, -57, -56, -55, -54, -53, -52, ...
                    -51.2, -51, -50, -45, -25.6, 0, 25.6] / 100;
%galvanic_coherences = [-100, -51.2, -30, -25.6, -20, -10, 0]/100;
galvanic_coherences = [-100, -65, -60, -59, -58, -57, -56, -55, -54, -53, ...
                       -52, -51.2, -51, -50, -45, -25.6, 0, 25.6] / 100;
anodic_coherences = fliplr([100, 65, 55, 51.2, 45, 40, 35, 34, 33, 32, 31, 30, ...
                            29, 28, 27, 26, 25.6, 12.8, 0]) / 100;
%pulse_coherences = [-100, -51.2, -25.6, 0, 25.6, 51.2, 100]/100;
%}
%{
control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6, 51.2] / 100;
pulse_coherences = [-100, -65, -60, -59, -58, -57, -56, -55, -51.2, -45, -25.6, 0, 25.6] / 100;
galvanic_coherences = [-100, -65, -60, -59, -58, -57, -56, -55, -51.2, -45, -25.6, 0, 25.6] / 100;
anodic_coherences = fliplr([100, 65, 55, 51.2, 45, 40, 35, 32, 31, 30, 29, 28, 25.6, 12.8, 0]) / 100;
%}
control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6, 51.2, 100] / 100;
pulse_coherences = [-100, -82.6, -69.8, -63.4, -60.2, -57, -53.8, -50.6, ...
                    -44.2, -31.4, -5.8, 0, 43, 100] ./ 100;
galvanic_coherences = [-100, -82.6, -69.8, -63.4, -60.2, -57, -53.8, -50.6, ...
                       -44.2, -31.4, -5.8, 0, 43, 100] ./ 100;
anodic_coherences = fliplr([100, 81.2, 55.6, 42.8, 36.4, 33.2, 30, 26.8, ...
                             23.6, 17.2, 4.4, 0, -21.2, -70, -100]) ./ 100;
pulse_amps = [-10*1e-6];
dc_amps = [-1.4, 0, 1.4]*1e-6;
%dc_amps = [-0.6, 0]*1e-6;
stim_amps = [pulse_amps, dc_amps];

%% Plot Firing Rates
ex_c = 0/100;
ex_trial = 19;
ex_stim_j = 3;
plot_name = "p1_only"; % 'single_stim' or 'subplot' or 'p1_only'
ylims = [];
plot_frs(sim_name, pulse_amps, stim_amps, p, f, N, N_E, t, t_task,...
                  t_taskoff, default_colors, ex_stim_j, ex_c, ex_trial, ylims, plot_name);


%% Plot Rasters
top_N = num_affected;%num_group*0.1;
ex_neurons = [3];
tlim = [0, 4];
ex_stim_j = 1;
ex_brain = 1;
ex_trial = 7;
ex_c = [-55, -55, 0, 35]/100;
plot_name = "subplot"; %'subplot' or 'ex_trial' or 'single_stim' or 'grouped_stim'
plot_rasters(sim_name, pulse_amps, stim_amps, ex_neurons, t, ...
                      t_task, t_taskoff, stim_freq, default_colors, top_N, ...
                      ex_stim_j, ex_c, ex_trial, tlim, plot_name);
                  
%% Plot Coefficeint of Variation (CV)
win_size = floor(0.250 / dt); %250ms moving window
cv_window = t >= 2.5 & t<3; %Plotting window
t_cut = 1; %omit trials with DTs longer than t_cut
ex_neuron = 7;
ex_trial = 1;
top_N = floor(1*num_group);
ex_c = [-55, -55, 0, 35]/100;
plot_name = "p1_wins"; %ex_neuron or 'ex_trial' or 'p1_wins'
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
plot_cv(sim_name, sim_names, pulse_amps, stim_amps, t, t_cut, N, top_N, num_group, ...
             num_affected, win_size, cv_window, default_colors, ex_c, ex_trial, ...
             ex_neuron, pulse_coherences, galvanic_coherences, control_coherences, ...
             anodic_coherences, start_trial, end_trial, num_trials, plot_name);


%% Plot Phaselocking to Pulses
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
idx_diff = stim_ind+1;% how far off timing is from pulse timing + 1 to account for t(1) = 0
win_start = 2.5;
win_stop = 3;
top_N = floor(1*num_group);
ex_c = [-55, -55, 0, 35] / 100;
t_cut = 1; %omit trials with DTs longer than t_cut
plot_phaselock(sim_names, pulse_amps, stim_amps, t, t_task, t_taskoff, t_cut, stim_freq, ...
                        num_group, num_affected, top_N, win_start, win_stop, idx_diff, ...
                        default_colors, start_trial, end_trial, num_trials, ex_c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences);

%% Plot Synchrony
N_start = 1;
N_end = floor(num_group);
win_start = 2.5;
win_stop = 3;
c_win = 300*1e-6;
ex_c = [-55, -55, 0, 35]/100;
sim_names = ["Brainless_m=0_Con", "Brainless_m=0_Discon"]; % [Connected, Disconnected]
symmetric = true;
t_cut = 1; %omit trials with DTs longer than t_cut
plot_sync(sim_names, pulse_amps, stim_amps, t, t_cut, num_group, num_affected, ...
                        N_start, N_end, win_start, win_stop, c_win, ex_c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        anodic_coherences, start_trial, end_trial, num_trials, symmetric);

%% Plot Firing Rate Distribution over distance from electrode
% stim_duration offset is used to account for onset spike of pulse trains
%win_start = t_task + stim_duration;
win_start = t_taskoff - 1/stim_freq + stim_duration - 0.1;

% omitting last 5ms of task period for exact comparison PR and FR
win_stop = t_taskoff - 1/stim_freq + stim_duration;
%win_stop = t_task + stim_duration + 0.1;

t_cut = 1.5; %omit trials with DTs longer than t_cut
%ex_c = [0, 0, 0, 0];
ex_c = [-57, -57, 0, 30]/100;
%  plot_name = 'ex_c' or  'ex_c_zoom' or 'p1_wins' or 'p1_loses'
plot_name = "p1_loses";
%sim_names = ["Brainless_m=0_Discon", "P1_Int", "P1_Rec"];
sim_names = ["iScience_Con"];
plot_frdist(sim_names, ex_c, pulse_amps, stim_amps, t, t_cut, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     anodic_coherences, start_trial, end_trial, num_trials, plot_name);

%% Plot Decisions and Decision Times
plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, ...
                num_batch, num_trials, pulse_coherences, ...
                galvanic_coherences, control_coherences, anodic_coherences);
            
%% Plot FR Trajectories
ex_c = [-55, -55, 0, 30] ./ 100;
t_cut = 1.5;
start_thresh = 10;
stop_thresh = 20;
plot_name = "p1_wins"; %p1_wins or p1_loses
plot_fr_trajectory(sim_name, pulse_amps, stim_amps, t, t_cut, t_task, ...
    ex_c, pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, N, p, f, N_E, ...
    start_thresh, stop_thresh, plot_name)

%% Plot FR Heatmaps
ex_c = [-55, -55, 0, 30] ./ 100;
t_cut = 2;
plot_name = "p1_loses"; %p1_wins or p1_loses
plot_fr_heatmap(sim_name, pulse_amps, stim_amps, t, t_cut, ex_c, ...
    pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, num_group, num_affected, plot_name)

%% Plot Pulsatile Blocking Validation
plot_name = "amp"; %"amp" or "spont"
ex_amp = -100e-6;
ex_spont = 1800;
plot_psval(sim_path, plot_name, pulse_amps, ex_amp, ex_spont)

%% Plot Galvanic Blocking Validation
save_amp = -1e-6;
plot_name = "external"; %"internal" or "external" (current amplitude for plotting)
plot_gsval(sim_path, save_amp, plot_name)
