%%% Paul Adkisson
%%% 9.6.21
%%% Purpose: Calculate decision time and accuracy from population firing
%%% rates
sim_name = "PR=500Hz_Discon";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))


num_batch = 3;
connected = false;

start_trial = 1;
end_trial = 28;
trials = start_trial:end_trial;
num_trials = length(trials);

control_coherences = 0;
galvanic_coherences = 0;
pulse_coherences = 0;
%control_coherences = [-100, -51.2, -25.6, 0, 25.6, 51.2] / 100;
%pulse_coherences = [-100, -65, -55, -51.2, -45, -25.6, 0, 25.6] / 100;
%galvanic_coherences = [-100, -65, -55, -51.2, -45, -25.6, 0, 25.6] / 100;
%galvanic_coherences = [100, 65, 55, 51.2, 45, 40, 35, 30, 25.6, 12.8, 0] / 100;

pulse_amps = [-50*1e-6];
%pulse_amps = [];
dc_amps = [-1.4, 0]*1e-6;
%dc_amps = [];
stim_amps = [pulse_amps, dc_amps];

for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
    if pulse
        fprintf("Pulse Stimulation Amplitude: %0.2fuA \n", stim_amp*1e6)
        output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
            [sim_name, stim_amp*1e6]);
        coherences = pulse_coherences;
    else
        fprintf("Galvanic Stimulation Amplitude: %0.2fuA \n", stim_amp*1e6)
        output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
            [sim_name, stim_amp*1e6]);
        if stim_amp == 0
            coherences = control_coherences;
        else
            coherences = galvanic_coherences;
        end
    end
    decisions = zeros(num_trials, length(coherences)); %0 for no decision, 1 for P1, 2 for P2
    decision_times = zeros(num_trials, length(coherences));
    avg_dts = zeros(length(coherences), 1);
    std_dts = zeros(length(coherences), 1);
    avg_acc = zeros(length(coherences), 1);
    batch_acc = zeros(length(coherences), num_batch);
    percent_nodec = zeros(length(coherences), 1);
    percent_earlydec = zeros(length(coherences), 1);
    for i = 1:length(coherences)
        c = coherences(i);
        fprintf("Coherence: %0.1f%% \n", c*100)
        output_coherentpath = strcat(output_stimpath, sprintf("/c=%0.3f", c));
        for trial = start_trial:end_trial
            output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));
            relative_trial = trial - start_trial + 1;
            load(output_trialpath, "recspikes")
            [pop_frs, fr_vars] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
            [decision, decision_idx] = get_decision(pop_frs);
            decisions(relative_trial, i) = decision;
            if decision_idx == -1
                decision_times(relative_trial, i) = NaN;
            else
                decision_times(relative_trial, i) = t(decision_idx) - t_task;
            end
            [dec, dec_idx] = get_decision(pop_frs);
            decisions(relative_trial, i) = dec;
            if dec_idx == -1
                decision_times(relative_trial, i) = NaN;
            else
                decision_times(relative_trial, i) = t(dec_idx) - t_task;
            end
        end
        %omit trials with no decision & those that decide before onset of task-related input
        coherent_times = decision_times(decisions(:, i)~=0 & decision_times(:, i)>0, i);
        avg_dts(i) = mean(coherent_times);
        std_dts(i) = std(coherent_times);
        %omit trials with no decision & those that decide before onset of task-related input
        coherent_decisions = decisions(decisions(:, i)~=0 & decision_times(:, i)>0, i);
        avg_acc(i) = sum(coherent_decisions == 1) / length(coherent_decisions);
        percent_nodec(i) = sum(decisions(:, i)==0, 'all') / end_trial;
        percent_earlydec(i) = sum(decision_times(:, i)<=0, 'all') / end_trial;

        %batch accuracy for statistical comparison
        batch_size = floor(length(coherent_decisions) / num_batch);
        for batch = 1:num_batch
            batch_idx = 1+(batch-1)*batch_size:batch*batch_size;
            batch_acc(i, batch) = sum(coherent_decisions(batch_idx)==1) / batch_size;
        end
    end
    if connected
        %Fit to logistic FN as in Hanks et al. 2006
        coeffs = lsqcurvefit(@logistic_acc, [1, 1], coherences, avg_acc');
        %Batch for statistical comparison
        batch_coeffs = zeros(num_batch, 2);
        for batch = 1:num_batch
            batch_coeffs(batch, :) = lsqcurvefit(@logistic_acc, [1, 1], coherences, batch_acc(:, batch)');
        end
    else
        coeffs = [];
        batch_coeffs = [];
    end
    decisionpath = strcat(output_stimpath, "/decisions.mat");
    save(decisionpath, "decisions", "decision_times", "avg_dts", "std_dts", ...
        "avg_acc", "percent_nodec", "coeffs", "batch_coeffs")
end

function [decision, dec_idx] = get_decision(pop_frs)
    decision_thresh = 15; %Hz
    if pop_frs(end, 1) < decision_thresh && pop_frs(end, 2) < decision_thresh || ...
            (pop_frs(end, 1) >= decision_thresh && pop_frs(end, 2) >= decision_thresh)%no decision
        decision = 0;
        dec_idx = -1;
    elseif pop_frs(end, 1) >= decision_thresh
        decision = 1;
        dec_idx = find(pop_frs(:, 1)>=decision_thresh, 1);
    else
        decision = 2;
        dec_idx = find(pop_frs(:, 2)>=decision_thresh, 1);
    end
end