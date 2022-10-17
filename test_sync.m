%%% Paul Adkisson
%%% 2.1.21
%%% Testing Synchrony Metric
c_win = 0.3*1e-3; %coincidence window (one-sided)
t = 0:0.05*1e-3:4;
sync_colors = [[0, 0, 0]; [1, 0, 0]];
%{
num_n1 = 1000;
num_n2 = 10;
jitter_factor = 2;
n1_spiketimes = rand(num_n1, 1)*4;
n2_spiketimes = n1_spiketimes(1:num_n2) + rand(num_n2, 1)*c_win*jitter_factor;
%}

load("Simulation EMBC Disconnected/bam_constants.mat")
%load("Simulation EMBC Disconnected/brain1/data/-10000.0nA_pulse/c=0.000/trial1.mat", "recspikes")
%load("Simulation EMBC Disconnected/brain1/data/0.0nA_galvanic/c=0.000/trial1.mat", "recspikes")
load("Simulation EMBC Disconnected/brain1/data/-28.0nA_galvanic/c=0.000/trial1.mat", "recspikes")
load("Simulation EMBC Disconnected/brain1/r.mat", "ball_r")
%}
%load("Simulation EMBC I_b100/bam_constants.mat")
%load("Simulation EMBC I_b100/brain1/data/-10000.0nA_pulse/c=0.000/trial4.mat", "recspikes")
%load("Simulation EMBC I_b100/brain1/data/-28.0nA_galvanic/c=0.000/trial4.mat", "recspikes")
%load("Simulation EMBC I_b100/brain1/data/0.0nA_galvanic/c=0.000/trial4.mat", "recspikes")
%load("Simulation EMBC I_b100/brain1/r.mat", "ball_r")
n1_spiketimes = t(recspikes("5"))';
n2_spiketimes = t(recspikes("6"))';
[sync, c_spikes] = get_sync(n1_spiketimes, n2_spiketimes, c_win);
figure;
hold on
scatter(n1_spiketimes, ones(length(n1_spiketimes), 1), 1000, sync_colors(c_spikes+1, :), 'Marker', '|')
scatter(n2_spiketimes, 2*ones(length(n2_spiketimes), 1), 1000, "k", 'Marker', '|')
ylim([0, 3])
hold off
title("Galvanic")

%{
%simulated rasters
avg_sync = 0;
num_iter = 1000;
for i = 1:num_iter
    n2_spiketimes = n1_spiketimes(1:num_n2) + rand(num_n2, 1)*c_win*jitter_factor;
    %[sync, ~] = get_sync(n1_spiketimes, n2_spiketimes, c_win);
    [sync, ~] = get_sync(n2_spiketimes, n1_spiketimes, c_win);
    avg_sync = avg_sync + sync/num_iter;
end
avg_sync
%}

N_start = 1;
N_end = floor(num_group*0.2);
win_start = 0.125;
win_stop = 0.375;
sync_win = 0.250; %250ms moving average window
sync_step = 0.010; %5ms moving average step size
t_down = sync_win:sync_step:t_span-sync_win;
pulse_timesync = zeros(length(t_down), 1);
galvanic_timesync = zeros(length(t_down), 1);
control_timesync = zeros(length(t_down), 1);
for stim = ["pulse", "galvanic", "control"]
    if stim == "pulse"
        load("Simulation EMBC I_b100/brain1/data/-10000.0nA_pulse/c=0.000/trial4.mat", "recspikes")
        pulse_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
    elseif stim == "galvanic"
        load("Simulation EMBC I_b100/brain1/data/-28.0nA_galvanic/c=0.000/trial4.mat", "recspikes")
        galvanic_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
    else
        load("Simulation EMBC I_b100/brain1/data/0.0nA_galvanic/c=0.000/trial4.mat", "recspikes")
        control_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
    end
    %{
    for i = 1:length(t_down)
        center_t = t_down(i)
        win_start = center_t - sync_win/2;
        win_stop = center_t + sync_win/2;
        sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
        if stim == "pulse"
            pulse_timesync(i) = mean(sync, 'all', 'omitnan');
        elseif stim == "galvanic"
            galvanic_timesync(i) = mean(sync, 'all', 'omitnan');
        else
            control_timesync(i) = mean(sync, 'all', 'omitnan');
        end
    end
    %}
end
figure;
hold on
plot(t_down, pulse_timesync, 'r')
plot(t_down, galvanic_timesync, 'g')
plot(t_down, control_timesync, 'k')
hold off

nan_color = uint8([0, 0, 128]);
figure;
h = heatmap(pulse_sync*100, 'ColorLimits', [0, 100], 'MissingDataColor', nan_color);
colormap(hot)
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
xlabel("neuron 2")
ylabel("neuron 1")
title("Pulsatile")

figure;
h = heatmap(galvanic_sync*100, 'ColorLimits', [0, 100], 'MissingDataColor', nan_color);
colormap(hot)
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
xlabel("neuron 2")
ylabel("neuron 1")
title("Galvanic")

figure;
h = heatmap(control_sync*100, 'ColorLimits', [0, 100], 'MissingDataColor', nan_color);
colormap(hot)
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
xlabel("neuron 2")
ylabel("neuron 1")
title("Control")

avg_pulse_sync = mean(pulse_sync, 'all', 'omitnan')*100
avg_galvanic_sync = mean(galvanic_sync, 'all', 'omitnan')*100
avg_control_sync = mean(control_sync, 'all', 'omitnan')*100