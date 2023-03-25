%%% Paul Adkisson
%%% 2/14/2023
%%% Plot Mean Firing Rate Heatmaps (over space) for each Stimulation Condition
function plot_fr_heatmap(sim_name, pulse_amps, stim_amps, t, t_cut, ex_c, ...
    pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, num_group, num_affected, plot_name)
    if plot_name == "p1_wins"
        win_num = 1;
    else
        win_num = 2;
    end
    dt = t(2) - t(1);
    stim_frs = zeros(length(stim_amps), num_trials, length(t), num_group);
    stim_cvs = zeros(length(stim_amps), num_trials, length(t));
    stim_skews = zeros(length(stim_amps), num_trials, length(t));
    stim_kurts = zeros(length(stim_amps), num_trials, length(t));
    load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r");
    ball_rs = get_ball_rs(ball_r, num_affected, num_group);
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        c = ex_c(j);
        pulse = j<=length(pulse_amps);
        if pulse
            disp("Pulsatile")
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
        [sim_name, stim_amp*1e6]);
        stim_coherences = pulse_coherences;
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
        [sim_name, stim_amp*1e6]);
            if stim_amp < 0 %cathodic GS
                disp("Cathodic GS")
                stim_coherences = galvanic_coherences;
            elseif stim_amp == 0 % control
                disp("Control")
                stim_coherences = control_coherences;
            else
                disp("Anodic GS")
                stim_coherences = anodic_coherences;
            end
        end
        load(strcat(output_stimpath, "/decisions.mat"), "decisions", "decision_times")
        for trial = start_trial:end_trial
            fprintf("Trial: %0.0f \n", trial)
            relative_trial = trial - start_trial + 1;
            dec_time = decision_times(relative_trial, stim_coherences==c);
            if decisions(relative_trial, stim_coherences==c) ~= win_num || ...
                    dec_time > t_cut
                stim_frs(j, relative_trial, :) = NaN;
                stim_cvs(j, relative_trial, :) = NaN;
                stim_skews(j, relative_trial, :) = NaN;
                stim_kurts(j, relative_trial, :) = NaN;
                continue %skip trials where P1 doesn't win/lose or decision takes too long
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                   "recspikes")
            win_size = 50e-3;
            avg_win_size = 200e-3;
            neuron_frs = recspikes2neuron_frs(recspikes, win_size, avg_win_size, t, num_group);
            stim_frs(j, trial, :, :) = neuron_frs;
            stim_cvs(j, trial, :) = std(neuron_frs, [], 2) ./ mean(neuron_frs, 2);
            stim_skews(j, trial, :) = skewness(neuron_frs, 0, 2); %flag = 0 to correct for sample bias
            stim_kurts(j, trial, :) = kurtosis(neuron_frs, 0, 2); %flag = 0 to correct for sample bias
        end
    end
    ps_frs = reshape(mean(stim_frs(1, :, :, :), 2, 'omitnan'), [length(t), num_group]);
    ps_cv = reshape(mean(stim_cvs(1, :, :), 2, 'omitnan'), [length(t), 1]);
    ps_skew = reshape(mean(stim_skews(1, :, :), 2, 'omitnan'), [length(t), 1]);
    ps_kurt = reshape(mean(stim_kurts(1, :, :), 2, 'omitnan'), [length(t), 1]);
    
    gs_frs = reshape(mean(stim_frs(2, :, :, :), 2, 'omitnan'), [length(t), num_group]);
    gs_cv = reshape(mean(stim_cvs(2, :, :), 2, 'omitnan'), [length(t), 1]);
    gs_skew = reshape(mean(stim_skews(2, :, :), 2, 'omitnan'), [length(t), 1]);
    gs_kurt = reshape(mean(stim_kurts(2, :, :), 2, 'omitnan'), [length(t), 1]);
    
    ctrl_frs = reshape(mean(stim_frs(3, :, :, :), 2, 'omitnan'), [length(t), num_group]);
    ctrl_cv = reshape(mean(stim_cvs(3, :, :), 2, 'omitnan'), [length(t), 1]);
    ctrl_skew = reshape(mean(stim_skews(3, :, :), 2, 'omitnan'), [length(t), 1]);
    ctrl_kurt = reshape(mean(stim_kurts(3, :, :), 2, 'omitnan'), [length(t), 1]);
    
    an_frs = reshape(mean(stim_frs(4, :, :, :), 2, 'omitnan'), [length(t), num_group]);
    an_cv = reshape(mean(stim_cvs(4, :, :), 2, 'omitnan'), [length(t), 1]);
    an_skew = reshape(mean(stim_skews(4, :, :), 2, 'omitnan'), [length(t), 1]);
    an_kurt = reshape(mean(stim_kurts(4, :, :), 2, 'omitnan'), [length(t), 1]);
    
    % Downsample for appropriate time-frequency resolution
    n_pts_per_bin = 4; %Number of time points per averaging window
    down_factor = floor(avg_win_size/(n_pts_per_bin*dt));
    t_down = t(1:down_factor:end);
    ball_rs = ball_rs(1:num_affected);
    ps_down = ps_frs(1:down_factor:end, 1:num_affected);
    gs_down = gs_frs(1:down_factor:end, 1:num_affected);
    ctrl_down = ctrl_frs(1:down_factor:end, 1:num_affected);
    an_down = an_frs(1:down_factor:end, 1:num_affected);
    
    xlims = [0.5, 3.5]; %omit first 0.5s to avoid 0FR artifact 
    figure;
    hold on
    plot(t, ps_cv, 'Color', default_colors(7, :))
    plot(t, gs_cv, 'Color', default_colors(5, :))
    plot(t, an_cv, 'Color', default_colors(6, :))
    plot(t, ctrl_cv, 'k-')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylabel("Coefficient of Variation (unitless)")
    hold off
    
    figure;
    hold on
    plot(t, ps_skew, 'Color', default_colors(7, :))
    plot(t, gs_skew, 'Color', default_colors(5, :))
    plot(t, an_skew, 'Color', default_colors(6, :))
    plot(t, ctrl_skew, 'k')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylabel("Skewness (unitless)")
    hold off
    
    figure;
    hold on
    plot(t, ps_kurt, 'Color', default_colors(7, : ))
    plot(t, gs_kurt, 'Color', default_colors(5, :))
    plot(t, an_kurt, 'Color', default_colors(6, :))
    plot(t, ctrl_kurt, 'k')
    legend(["PS", "CGS", "AGS", "No Stim"])
    xlabel("Time (s)")
    xlim(xlims)
    ylabel("Kurtosis (unitless)")
    hold off
    
    near_mask = ball_rs<=500e-6;
    %far_mask = ~near_mask;
    far_mask = logical(ones(size(ball_rs)));
    near_ticks = 0:100:500;
    %far_ticks = 500:500:2000;
    far_ticks = 0:500:2000;
    near_map = pink;
    far_map = hot;
    near_lims = [0, 150];
    far_lims = [0, 15];
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, ps_down(:, near_mask)');
    colormap(near_map)
    caxis([0, 150])
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("Pulse Near")
    
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ps_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Pulse Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, gs_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("CGS Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, gs_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("CGS Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, ctrl_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("Control Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, ctrl_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("Control Far")
    
    figure;
    s = pcolor(t_down, ball_rs(near_mask)*1e6, an_down(:, near_mask)');
    colormap(near_map)
    caxis(near_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(near_ticks)
    title("AGS Near")
    
    figure;
    s = pcolor(t_down, ball_rs(far_mask)*1e6, an_down(:, far_mask)');
    colormap(far_map)
    caxis(far_lims)
    cbar = colorbar('FontSize', 20);
    cbar.Label.String = "P1 Firing Rates (spk/s)";
    set(s, 'EdgeColor', 'none');
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    yticks(far_ticks)
    title("AGS Far")
end

function neuron_frs = recspikes2neuron_frs(recspikes, win_size, avg_win_size, t, N)
    dt = t(2) - t(1);
    spikes = zeros(length(t), N);
    for nn = 1:N
        for spikeidx = recspikes(int2str(nn))
            spikes(spikeidx, nn) = 1;
        end
    end
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
end