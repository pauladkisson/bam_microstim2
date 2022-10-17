clear;
sim_num = 2;
sim_path = sprintf("Simulation %0.0f", sim_num);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");

%{
%Example FRs
ex_c = -51.2/100;
ex_trial = 1;
%ex_stim_amp = -10*1e-6;
ex_stim_amp = -1*1e-9;
pulse = false;
if pulse
    output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, ex_stim_amp*1e9]);
else
    output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
        [sim_num, ex_stim_amp*1e9]);
end
load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "pop_frs")
figure;
hold on
for i = 1:p+2
    plot(t, pop_frs(:, i))
end
hold off
legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
title(sprintf("Coherence: %0.1f%%, Trial %0.0f", [ex_c*100, ex_trial]))
%}

%{
%Example Raster by pop-type and total activity 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = -1;
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if sim_num >= 2
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, stim_amp*1e9]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
                [sim_num, stim_amp*1e9]);
        end
    else
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fuA_pulse", [sim_num, stim_amp*1e6]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fpA_galvanic", ...
                [sim_num, stim_amp*1e12]);
        end
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "recspikes", "pop_frs")
    spikes = zeros(length(t), N);
    neuron_num = 1:N;
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    
    g1_spikes = spikes(:, 1:num_group);
    [g1_time_idx, g1_neuron_idx] = get_spike_idx(g1_spikes);
    g2_spikes = spikes(:, num_group+1:2*num_group);
    [g2_time_idx, g2_neuron_idx] = get_spike_idx(g2_spikes);
    ns_spikes = spikes(:, 2*num_group+1:N_E);
    [ns_time_idx, ns_neuron_idx] = get_spike_idx(ns_spikes);
    int_spikes = spikes(:, N_E+1:end);
    [int_time_idx, int_neuron_idx] = get_spike_idx(int_spikes);
    
    figure;
    hold on
    scatter(t(g1_time_idx), g1_neuron_idx, "Marker", "|", ...
        "MarkerFaceColor", default_colors(1, :), "MarkerEdgeColor", default_colors(1, :))
    scatter(t(g2_time_idx), num_group+g2_neuron_idx, "Marker", "|", ...
        "MarkerFaceColor", default_colors(2, :), "MarkerEdgeColor", default_colors(2, :))
    scatter(t(ns_time_idx), 2*num_group+ns_neuron_idx, "Marker", "|", ...
        "MarkerFaceColor", default_colors(3, :), "MarkerEdgeColor", default_colors(3, :))
    scatter(t(int_time_idx), N_E+int_neuron_idx, "Marker", "|", ...
        "MarkerFaceColor", default_colors(4, :), "MarkerEdgeColor", default_colors(4, :))
    xlabel("Time (s)")
    ylabel("Ranked Neuron Activity")
    legend(["Population 1", "Population 2", "Non-Selective", "Interneurons"])
    if sim_num >= 2
        if pulse
            title(sprintf("%0.0fnA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        else
            title(sprintf("%0.0fnA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        end
    else
        if pulse
            title(sprintf("%0.0fuA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e6, ex_c*100]))
        else
            title(sprintf("%0.0fpA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e12, ex_c*100]))
        end
    end
end
%}

%{
%Example FRs with variance 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = -0.512;
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if sim_num >= 2
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, stim_amp*1e9]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
                [sim_num, stim_amp*1e9]);
        end
    else
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fuA_pulse", [sim_num, stim_amp*1e6]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fpA_galvanic", ...
                [sim_num, stim_amp*1e12]);
        end
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "recspikes", "pop_frs")
    spikes = zeros(length(t), N);
    neuron_num = 1:N;
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    [pop_frs, fr_vars] = spikes2popfrs(spikes, dt, p, f, N_E);
    figure;
    hold on
    n_down = 20;
    for i = 1:4
        shadedErrorBar(t(1:n_down:end) , pop_frs(1:n_down:end, i), fr_vars(1:n_down:end, i), ...
            'lineProps', {'Color', default_colors(i, :)})
    end
    hold off
    xlabel("Time (s)")
    ylabel("Firing Rate (Hz)")
    legend(["Population 1", "Population 2", "Non-Selective", "Interneurons"])
    if sim_num >= 2
        if pulse
            title(sprintf("%0.0fnA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        else
            title(sprintf("%0.0fnA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        end
    else
        if pulse
            title(sprintf("%0.0fuA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e6, ex_c*100]))
        else
            title(sprintf("%0.0fpA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e12, ex_c*100]))
        end
    end
end
%}

%Total Spike Histograms
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = 1;
ex_c_idx = coherences == ex_c;
figure;
hold on
axs = [];
binEdges = 0:0.5:60;
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if sim_num >= 2
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, stim_amp*1e9]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
                [sim_num, stim_amp*1e9]);
        end
    else
        if pulse
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fuA_pulse", [sim_num, stim_amp*1e6]);
        else
            output_stimpath = sprintf("Simulation %0.0f/data/%0.0fpA_galvanic", ...
                [sim_num, stim_amp*1e12]);
        end
    end
    load(strcat(output_stimpath, "/decisions.mat"), ...
        "totspikes_g1", "totspikes_g2", "totspikes_ns", "totspikes_int")
    axs(j) = subplot(3, 1, j);
    histogram(totspikes_g1(ex_c_idx, :), 'BinEdges', binEdges);
    if pulse
        title(sprintf("%0.0fnA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
    else
        title(sprintf("%0.0fnA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
    end
    if j == length(stim_amps)
        xlabel("Average Firing Rate (Hz)")
    elseif j == 2
        ylabel("Number of Neurons")
    end
end
linkaxes(axs);
hold off
%}


%{
%Decision Times
figure;
hAx = axes;
%hAx.XScale = 'log';
hold on
fig_leg = [];
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
    if pulse
        datapath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", ...
            [sim_num, stim_amp*1e9]);
        fig_leg = [fig_leg, sprintf("%0.0fnA Pulse", stim_amp*1e9)];
    else
        datapath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
            [sim_num, stim_amp*1e9]);
        fig_leg = [fig_leg, sprintf("%0.0fnA Galvanic", stim_amp*1e9)];
    end
    load(strcat(datapath, "/decisions.mat"), "avg_dts", "std_dts");
    errorbar(coherences, avg_dts, std_dts)
end
hold off
xlabel("Coherence")
ylabel("Decision Time (s)")
legend(fig_leg)
%}

%{
%Accuracies
c = 0:0.01:1;
figure;
hAx = axes;
%hAx.XScale = 'log';
hold on
idx = 1;
fig_leg = [];
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
    if pulse
        datapath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", ...
            [sim_num, stim_amp*1e9]);
        fig_leg = [fig_leg, sprintf("%0.0fnA Pulse", stim_amp*1e9)];
    else
        datapath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
            [sim_num, stim_amp*1e9]);
        fig_leg = [fig_leg, sprintf("%0.0fnA Galvanic", stim_amp*1e9)];
    end
    load(strcat(datapath, "/decisions.mat"), "avg_acc", "decisions");
    %coeffs
    %scatter(coherences, avg_acc, dc_colors(idx))
    %plot(c, weibull(coeffs, c), dc_colors(idx))
    plot(coherences, avg_acc, 'o-')
end
hold off
xlabel("Coherence")
ylabel("Accuracy")
legend(fig_leg)
%f = flipud(get(gca, 'Children'));
%legend([f(2), f(4), f(6)], "I-dc=-4pA", "I-dc=0pA", "I-dc=4pA")
%legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%{
%example trial schreiber correlation
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
corr_window = t>=2.5 & t<3;
ex_neuron = 2;
w = gausswin(floor(0.001/dt));
cross_corr = zeros(length(stim_amps), num_group);
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if pulse
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
            [sim_name, ex_brain, stim_amp*1e9]);
        stim_color = default_colors(7, :);
    else
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
            [sim_name, ex_brain, stim_amp*1e9]);
        if stim_amp == 0
            stim_color = [0, 0, 0];
        else
            stim_color = default_colors(5, :);
        end
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "recspikes")
    spikes = zeros(length(t), num_group);
    for nn = 1:num_group
        spikeidx = recspikes(int2str(nn));
        spikes(spikeidx, nn) = 1;
    end
    win_spikes = spikes(corr_window, :);
    S = zeros(sum(corr_window)+length(w)-1, num_group);
    for nn = 1:num_group
        S(:, nn) = conv(win_spikes(:, nn), w);
    end
    size(sum(S, 1)' * sum(S, 1))
    R = S' * S ./ (sum(S, 1)' * sum(S, 1));
    cross_corr(j, :) = mean(R, 1, 'omitnan');
end
figure;
hold on
scatter(ball_r*1e6, cross_corr(1, :), [], default_colors(7, :), 'filled')
scatter(ball_r*1e6, cross_corr(2, :), [], default_colors(5, :), 'filled')
scatter(ball_r*1e6, cross_corr(3, :), [], "k", 'filled')
hold off
ylabel("Schreiber Similarity")
xlabel("Distance from electrode (um)")
legend(["Pulsatile", "Galvanic", "Control"])
title("Disconnected")
%}

%Spike Timing Correlation 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = 0.256;
ex_brain = 1;
timeidx = t >= t_task & t < t_taskoff;
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if pulse
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    else
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "recspikes")
    spikes = zeros(length(t), N);
    neuron_num = 1:N;
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    fprintf("Stimulation %0.1fnA \n", stim_amp*1e9)
    neuron_frs = spikes2neuron_frs(spikes, dt);
    
    g1_frs = neuron_frs(timeidx, 1:num_group);
    g1_frs = g1_frs(:, ~all(g1_frs==0, 1)); %omit neurons that never fire
    g1_corr = cov(g1_frs) ./ sqrt(var(g1_frs) .* var(g1_frs)');
    g1_corr = g1_corr(~eye(size(g1_corr))); %omit diagonal
    g1_avg_corr = mean(g1_corr);
    fprintf("g1 : %0.4f \n", g1_avg_corr)
    
    g2_frs = neuron_frs(timeidx, num_group+1:2*num_group);
    g2_frs = g2_frs(:, ~all(g2_frs==0, 1)); %omit neurons that never fire
    g2_corr = cov(g2_frs) ./ sqrt(var(g2_frs) .* var(g2_frs)');
    g2_corr = g2_corr(~eye(size(g2_corr))); %omit diagonal
    g2_avg_corr = mean(g2_corr);
    fprintf("g2 : %0.4f \n", g2_avg_corr)
    
    ns_frs = neuron_frs(timeidx, 2*num_group+1:N_E);
    ns_frs = ns_frs(:, ~all(ns_frs==0, 1)); %omit neurons that never fire
    ns_corr = cov(ns_frs) ./ sqrt(var(ns_frs) .* var(ns_frs)');
    ns_corr = ns_corr(~eye(size(ns_corr))); %omit diagonal
    ns_avg_corr = mean(ns_corr);
    fprintf("ns : %0.4f \n", ns_avg_corr)
    
    int_frs = neuron_frs(timeidx, N_E+1:end);
    int_frs = int_frs(:, ~all(int_frs==0, 1)); %omit neurons that never fire
    int_corr = cov(int_frs) ./ sqrt(var(int_frs) .* var(int_frs)');
    int_corr = int_corr(~eye(size(int_corr))); %omit diagonal
    int_avg_corr = mean(int_corr);
    fprintf("int : %0.4f \n", int_avg_corr)
    
end

function [time_idx, neuron_idx] = get_spike_idx(g_spikes)
    [~, g_idx] = sort(sum(g_spikes, 1));
    g_spikes = g_spikes(:, g_idx);
    [time_idx, neuron_idx] = find(g_spikes);
end

function [frs, fr_var] = spikes2popfrs(spikes, dt, p, f, N_E)
    win_size = 5e-3;
    avg_win_size = 50e-3;
    num_group = floor(f*N_E);
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
    frs = zeros(size(spikes, 1), p+2);
    fr_var = zeros(size(spikes, 1), p+2);
    for i = 1:p
        group_idx = ((i-1)*num_group+1):i*num_group;
        frs(:, i) = mean(neuron_frs(:, group_idx), 2);
        fr_var(:, i) = std(neuron_frs(:, group_idx), [], 2);
    end
    frs(:, end-1) = mean(neuron_frs(:, f*p*N_E+1:N_E), 2);
    fr_var(:, end-1) = std(neuron_frs(:, f*p*N_E+1:N_E), [], 2);
    frs(:, end) = mean(neuron_frs(:, N_E+1:end), 2);
    fr_var(:, end) = std(neuron_frs(:, N_E+1:end), [], 2);
end