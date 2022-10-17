%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Synchrony
function plot_sync(sim_names, pulse_amps, stim_amps, t, num_group, ...
                        brains, num_brains, N_start, N_end, ...
                        win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials, default_colors)
    ps_norm_means = zeros(num_brains, length(sim_names));
    gs_norm_means = zeros(num_brains, length(sim_names));
    percent_phaselocked = 26.8/100;
    num_phaselocked = floor(percent_phaselocked * num_group);
    for sim_name = sim_names
        disp(sim_name)
        
        ball_rs = zeros(num_brains, num_group);
        stim_sync = zeros(length(stim_amps), num_brains, num_group, num_group);
        max_ps_sync = zeros(num_brains, num_group);
        mean_ps_sync = zeros(num_brains, 1);
        max_gs_sync = zeros(num_brains, num_group);
        mean_gs_sync = zeros(num_brains, 1);
        ps_phaselocked_sync = zeros(num_brains, 1);
        ps_unlocked_sync = zeros(num_brains, 1);
        for brain = brains
            fprintf("brain %0.0f \n", brain)
            load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, brain]), "ball_r")
            ball_rs(brain, :) = ball_r;
            for j = 1:length(stim_amps)
                stim_amp = stim_amps(j);
                pulse = j<=length(pulse_amps);
                if pulse
                    output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                        [sim_name, brain, stim_amp*1e9]);
                    stim_coherences = pulse_coherences;
                    disp("Pulsatile")
                else
                    output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                        [sim_name, brain, stim_amp*1e9]);
                    if stim_amp == 0
                        stim_coherences = control_coherences;
                        ctrl_ball_r = ball_r;
                        if brain ~= 1
                            continue
                        end
                        disp("Control")
                    else
                        stim_coherences = galvanic_coherences;
                        disp("Galvanic")
                    end
                end
                %load(strcat(output_stimpath, "/decisions.mat"), "final_decisions")
                load(sprintf("Simulation %s/brain1/data/0.0nA_galvanic/decisions.mat", ...
                        [sim_name]), 'final_decisions'); %using only control final decs
                if sim_name == sim_names(1)
                    num_wins = sum(final_decisions(start_trial:end_trial, stim_coherences==0)==1, 'all') * ones(num_group, num_group);
                else
                    num_wins = num_trials * ones(num_group, num_group);
                end
                nan_neurons = zeros(num_trials, num_group);
                for trial = start_trial:end_trial
                    fprintf("Trial %0.0f \n", trial)
                    if sim_name==sim_names(1) && final_decisions(trial, stim_coherences==c) ~= 1
                        disp("P1 lost")
                        continue %skip trials where P1 doesn't win
                    end
                    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                        "recspikes")
                    %pairwise_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                    pairwise_sync = get_sym_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                    no_spike_neurons = all(isnan(pairwise_sync), 1);
                    num_wins(no_spike_neurons, :) = num_wins(no_spike_neurons, :) - 1;
                    num_wins(:, no_spike_neurons) = num_wins(:, no_spike_neurons) - 1;
                    pairwise_sync = pairwise_sync(~no_spike_neurons, ~no_spike_neurons);
                    stim_sync(j, brain, ~no_spike_neurons, ~no_spike_neurons) = ...
                        reshape(stim_sync(j, brain, ~no_spike_neurons, ~no_spike_neurons), size(pairwise_sync)) + pairwise_sync;
                    nan_neurons(trial, :) = no_spike_neurons;
                end
                if any(all(nan_neurons==1, 1)) %no spikes for all 36 trials
                    stim_sync(j, brain, all(nan_neurons==1, 1), :) = NaN;
                    stim_sync(j, brain, :, all(nan_neurons==1, 1)) = NaN;
                end
                stim_sync(j, brain, :, :) = reshape(stim_sync(j, brain, :, :), size(num_wins)) ./ num_wins;
            end
            max_ps_sync(brain, :) = reshape(max(stim_sync(1, brain, :, :), [], 4), [1, num_group]);
            max_gs_sync(brain, :) = reshape(max(stim_sync(2, brain, :, :), [], 4), [1, num_group]);
            mean_ps_sync(brain) = reshape(mean(stim_sync(1, brain, :, :), 'all', 'omitnan'), [1, 1]);
            mean_gs_sync(brain) = reshape(mean(stim_sync(2, brain, :, :), 'all', 'omitnan'), [1, 1]);
            ps_phaselocked_sync(brain) = reshape(mean(stim_sync(1, brain, 3:num_phaselocked, 3:num_phaselocked), 'all', 'omitnan'), [1, 1]);
            ps_unlocked_sync(brain) = reshape(mean(stim_sync(1, brain, num_phaselocked+1:end, num_phaselocked+1:end), 'all', 'omitnan'), [1, 1]);
        end

        pulse_sync = reshape(mean(stim_sync(1, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        galvanic_sync = reshape(mean(stim_sync(2, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        control_sync = reshape(stim_sync(3, 1, :, :), [num_group, num_group]);

        if sim_name == "EMBC I_b100"
            save("matdata/connected.mat", "pulse_sync", "galvanic_sync", "control_sync")
        else
            save("matdata/disconnected.mat", "pulse_sync", "galvanic_sync", "control_sync")
        end

        nan_color = uint8([0, 0, 128]);
        ticks = log10([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]);
        tick_labels = ["1", "10", "", "", "", "", "", "", "", "", "100"];
        figure;
        h = heatmap(log10(pulse_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Pulsatile")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(galvanic_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Galvanic")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(control_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Control")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        
        %Calculate fraction of neuron pairs that are "significantly
        %synchronized" (>control max)
        sync_num = 1;
        ctrl_max = max(control_sync, [], 'all')
        ctrl_max = 0.037;
        ctrl_synced = control_sync > 0.037; %control synchronized above disconnected ctrl_max
        ctrl_num_synced = sum(ctrl_synced, 'all') / 2
        
        ps_perc_synced = zeros(num_brains, 1);
        ps_perc_atleastX = zeros(num_brains, 1);
        ps_perc_phase_n_sync = zeros(num_brains, 1);
        ps_maxsync = zeros(num_brains, 1);
        ps_num_synced = zeros(num_brains, 1);
        ps_d_synced = zeros(num_brains, 1);
        ps_maxsync_of_synced = zeros(num_brains, 1);
        ps_minsync_of_synced = zeros(num_brains, 1);
        ps_num_sync_of_synced = zeros(num_brains, 1);
        
        gs_perc_synced = zeros(num_brains, 1);
        gs_perc_atleastX = zeros(num_brains, 1);
        gs_maxsync = zeros(num_brains, 1);
        gs_num_synced = zeros(num_brains, 1);
        gs_d_synced = zeros(num_brains, 1);
        gs_maxsync_of_synced = zeros(num_brains, 1);
        gs_minsync_of_synced = zeros(num_brains, 1);
        gs_num_sync_of_synced = zeros(num_brains, 1);
        for brain = brains
            ps_synced = stim_sync(1, brain, :, :) > ctrl_max;
            ps_maxsync_of_synced(brain) = max(stim_sync(1, brain, ps_synced), [], 'all');
            ps_minsync_of_synced(brain) = min(stim_sync(1, brain, ps_synced), [], 'all');
            ps_perc_phase_n_sync(brain) = sum(ps_synced(:, :, 1:num_phaselocked, 1:num_phaselocked), 'all') / sum(ps_synced, 'all');
            ps_atleastX_synced = sum(ps_synced, 4) >= sync_num;
            ps_perc_synced(brain) = sum(ps_synced, 'all') / (length(ps_synced)^2 - num_group);
            ps_num_synced(brain) = sum(ps_synced, 'all') / 2;
            ps_perc_atleastX(brain) = sum(ps_atleastX_synced, 'all') / num_group;
            ps_maxsync(brain) = max(stim_sync(1, brain, :, :), [], 'all');
            ps_d_synced(brain) = max(ball_r(ps_atleastX_synced));
            ps_num_sync_of_synced(brain) = max(sum(ps_synced, 4), [], 'all');
            
            gs_synced = stim_sync(2, brain, :, :) > ctrl_max;
            if ~all(gs_synced==0, 'all')
                gs_maxsync_of_synced(brain) = max(stim_sync(2, brain, gs_synced), [], 'all');
                gs_minsync_of_synced(brain) = min(stim_sync(2, brain, gs_synced), [], 'all');
            end
            gs_atleastX_synced = sum(gs_synced, 4) >= sync_num;
            gs_perc_synced(brain) = sum(gs_synced, 'all') / (length(gs_synced)^2 - num_group);
            gs_num_synced(brain) = sum(gs_synced, 'all') / 2;
            gs_perc_atleastX(brain) = sum(gs_atleastX_synced, 'all') / num_group;
            gs_maxsync(brain) = max(stim_sync(2, brain, :, :), [], 'all');
            gs_num_sync_of_synced(brain) = max(sum(gs_synced, 4), [], 'all');
            try
                gs_d_synced(brain) = max(ball_r(gs_atleastX_synced));
            catch
                assert(all(gs_atleastX_synced==0, 'all')) % NaN if no synchronized connections
                gs_d_synced(brain) = NaN;
            end
        end
        [~, p_perc_synced] = ttest2(ps_perc_synced, gs_perc_synced)
        [~, p_perc_atleastX] = ttest2(ps_perc_atleastX, gs_perc_atleastX)
        
        mean_ps_sync = mean(ps_perc_synced)
        sem_ps_sync = std(ps_perc_synced) / sqrt(num_brains)
        mean_ps_atleastX = mean(ps_perc_atleastX)
        sem_ps_atleastX = std(ps_perc_atleastX) / sqrt(num_brains)
        mean_ps_perc_phase_n_sync = mean(ps_perc_phase_n_sync)
        sem_ps_perc_phase_n_sync = std(ps_perc_phase_n_sync) / sqrt(num_brains)
        mean_ps_maxsync = mean(ps_maxsync)
        sem_ps_maxsync = std(ps_maxsync) / sqrt(num_brains)
        mean_ps_num_synced = mean(ps_num_synced)
        sem_ps_num_synced = std(ps_num_synced) / sqrt(num_brains)
        max_ps_d_synced = max(ps_d_synced) * 1e6
        max_ps_sync_o_sync = max(ps_maxsync_of_synced)
        min_ps_sync_o_sync = min(ps_minsync_of_synced)
        max_ps_num_sync_o_sync = max(ps_num_sync_of_synced)
        
        mean_norm_gs_sync = mean(gs_perc_synced)
        sem_norm_gs_sync = std(gs_perc_synced) / sqrt(num_brains)
        mean_norm_gs_atleastX = mean(gs_perc_atleastX)
        sem_norm_gs_atleastX = std(gs_perc_atleastX) / sqrt(num_brains)
        mean_gs_maxsync = mean(gs_maxsync)
        sem_gs_maxsync = std(gs_maxsync) / sqrt(num_brains)
        mean_gs_num_synced = mean(gs_num_synced)
        sem_gs_num_synced = std(gs_num_synced) / sqrt(num_brains)
        max_gs_d_synced = max(gs_d_synced) * 1e6
        max_gs_sync_o_sync = max(gs_maxsync_of_synced)
        min_gs_sync_o_sync = min(gs_minsync_of_synced)
        max_gs_num_sync_o_sync = max(gs_num_sync_of_synced)
        
        %save for comparison connected vs disconnected
        if contains(sim_name, "Discon")
            discon_perc_synced = ps_perc_synced;
            discon_perc_atleastX = ps_perc_atleastX;
            discon_perc_syncNphase = ps_perc_phase_n_sync;
        else
            con_perc_synced = ps_perc_synced;
            con_perc_atleastX = ps_perc_atleastX;
            con_perc_syncNphase = ps_perc_phase_n_sync;
        end

        %average for plotting
        ps_synced = double(pulse_sync > ctrl_max);
        gs_synced = double(galvanic_sync > ctrl_max);
        figure;
        h = heatmap(ps_synced);
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        colormap(hot)
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Pulse")

        figure;
        h = heatmap(gs_synced);
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        colormap(hot)
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Galvanic")
        
        total_N = num_brains*num_group;
        tot_ps_sync = reshape(max_ps_sync, [total_N, 1]);
        tot_gs_sync = reshape(max_gs_sync, [total_N, 1]);
        max_ctrl_sync = max(control_sync, [], 2);
        tot_ball_rs = reshape(ball_rs, [total_N, 1]);
        if contains(sim_name, "Discon")
            shape = "o";
        else
            shape = "^";
        end
        figure(123);
        set(gca, 'fontsize', 18)
        hold on
        scatter(tot_ball_rs*1e6, tot_ps_sync*100, [], ones(total_N, 3).*default_colors(7, :), 'filled', shape)
        scatter(tot_ball_rs*1e6, tot_gs_sync*100, [], ones(total_N, 3).*default_colors(5, :), 'filled', shape)
        scatter(ctrl_ball_r*1e6, max_ctrl_sync*100, [], "k", 'filled', shape)
        hold off
        xlabel("Distance from Electrode (um)")
        ylabel("Maximum Spike Synchrony (%)")
        
        %Statistics
        %ctrl_mean = mean(max_ctrl_sync, 'omitnan')
        %ctrl_mean = mean(control_sync, 'all', 'omitnan')
        %ctrl_std = std(max_ctrl_sync, 'omitnan');
        %to allow comparison connected/disconnected, always normalize to
        %disconnected control
        ctrl_mean = 0.0139;
        %pulse_norm_mean = mean(max_ps_sync, 2) - ctrl_mean;
        %galvanic_norm_mean = mean(max_gs_sync, 2) - ctrl_mean;
        
        pulse_norm_mean = mean_ps_sync - ctrl_mean;
        phaselocked_norm_mean = ps_phaselocked_sync - ctrl_mean;
        unlocked_norm_mean = ps_unlocked_sync - ctrl_mean;
        galvanic_norm_mean = mean_gs_sync - ctrl_mean;
        stim_means = [mean(galvanic_norm_mean); mean(pulse_norm_mean)];
        stim_sems = [std(galvanic_norm_mean)/sqrt(num_brains); std(pulse_norm_mean)/sqrt(num_brains)];

        %[~, p] = ttest2(galvanic_norm_mean, pulse_norm_mean)
        [~, ~, stats] = anova1([phaselocked_norm_mean, unlocked_norm_mean, galvanic_norm_mean]);
        p = multcompare(stats)
        [~, p_gs] = ttest(galvanic_norm_mean)
        [~, p_unlocked] = ttest(unlocked_norm_mean)
        
        galvanic_norm_bar = mean(galvanic_norm_mean)
        galvanic_norm_sem = std(galvanic_norm_mean) ./ sqrt(length(galvanic_norm_mean))
        pulse_norm_bar = mean(pulse_norm_mean)
        pulse_norm_sem = std(pulse_norm_mean) ./ sqrt(length(pulse_norm_mean))
        phaselocked_norm_bar = mean(phaselocked_norm_mean)
        phaselocked_norm_sem = std(phaselocked_norm_mean) ./ sqrt(length(phaselocked_norm_mean))
        unlocked_norm_bar = mean(unlocked_norm_mean)
        unlocked_norm_sem = std(unlocked_norm_mean) ./ sqrt(length(unlocked_norm_mean))
        
        
        ps_norm_means(:, sim_name==sim_names) = pulse_norm_mean;
        gs_norm_means(:, sim_name==sim_names) = galvanic_norm_mean;
    end
    %More stats
    %[~, ~, stats] = anova1([ps_norm_means, gs_norm_means]);
    [~, ~, stats] = anova2([ps_norm_means; gs_norm_means], num_brains);
    c = multcompare(stats)
    c = multcompare(stats, 'Estimate', 'row')
    ps_norm_means
    gs_norm_means
    
    %More stats
    [~, p_sync] = ttest2(con_perc_synced, discon_perc_synced)
    [~, p_atleastX] = ttest2(con_perc_atleastX, discon_perc_atleastX)
    [~, p_syncNphase] = ttest2(con_perc_syncNphase, discon_perc_syncNphase)
end