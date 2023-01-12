%%% Paul Adkisson
%%% 12/6/2022
%%% Plot Mean Firing Rate Trajectories for each Stimulation Condition
function plot_fr_trajectory(sim_name, pulse_amps, stim_amps, t, t_cut, ex_c, ...
    pulse_coherences, galvanic_coherences, control_coherences, anodic_coherences, ...
    default_colors, start_trial, end_trial, num_trials, N, p, f, N_E, ...
    start_thresh, stop_thresh, plot_name)
    if plot_name == "p1_wins"
        win_num = 1;
    else
        win_num = 2;
    end
    dt = t(2) - t(1);
    stim_frs = zeros(length(stim_amps), num_trials, length(t));
    stim_slopes = zeros(length(stim_amps), num_trials);
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
                stim_slopes(j, relative_trial) = NaN;
                continue %skip trials where P1 doesn't win/lose or decision takes too long
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                   "recspikes")
            [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
            stim_frs(j, relative_trial, :) = pop_frs(:, 1);
            if plot_name == "p1_wins"
                preidx = find(pop_frs(:, 1)>=start_thresh, 1);
                postidx = find(pop_frs(:, 1)>=stop_thresh, 1);
                slope_t = t(preidx:postidx);
                slope = (stop_thresh-start_thresh) / (t(postidx)-t(preidx));
                coeffs = [slope, pop_frs(preidx, 1) - slope*t(preidx)];
                stim_slopes(j, relative_trial) = coeffs(1);
                
                %plot to debug outliers
                debug = slope > 100; 
                if debug
                    figure;
                    slope_y = coeffs(2) + coeffs(1)*slope_t;
                    hold on;
                    plot(t, pop_frs(:, 1))
                    plot(slope_t, slope_y, "r--")
                    title(sprintf("j=%0.0f", j))
                    %break
                end                     
            end
        end
    end
    pulse_frs = reshape(stim_frs(1, :, :), [num_trials, length(t)]);
    galvanic_frs = reshape(stim_frs(2, :, :), [num_trials, length(t)]);
    control_frs = reshape(stim_frs(3, :, :), [num_trials, length(t)]);
    anodic_frs = reshape(stim_frs(4, :, :), [num_trials, length(t)]);
    pulse_trialmean = mean(pulse_frs, 1, 'omitnan');
    galvanic_trialmean = mean(galvanic_frs, 1, 'omitnan');
    control_trialmean = mean(control_frs, 1, 'omitnan');
    anodic_trialmean = mean(anodic_frs, 1, 'omitnan');
    
    figure;
    set(gca, 'Fontsize', 18)
    hold on
    plot(t, pulse_trialmean, 'Color', default_colors(7, :), 'Linewidth', 2)
    plot(t, galvanic_trialmean, 'Color', default_colors(5, :), 'Linewidth', 2)
    plot(t, control_trialmean, "k", 'Linewidth', 2)
    plot(t, anodic_trialmean, 'Color', default_colors(6, :), 'Linewidth', 2)
    hold off
    xlabel("Time (s)")
    ylabel("P1 Firing Rate (spk/s)")
    if plot_name == "p1_wins"
        title("P1 Wins")
    else
        title("P1 Loses")
    end
    
    pulse_slopes = reshape(stim_slopes(1, :, :), [num_trials, 1]);
    galvanic_slopes = reshape(stim_slopes(2, :, :), [num_trials, 1]);
    control_slopes = reshape(stim_slopes(3, :, :), [num_trials, 1]);
    anodic_slopes = reshape(stim_slopes(4, :, :), [num_trials, 1]);
    mean_ctrl = mean(control_slopes, 'omitnan');
    norm_ps = pulse_slopes - mean_ctrl;
    norm_gs = galvanic_slopes - mean_ctrl;
    norm_an = anodic_slopes - mean_ctrl;
    mean_slopes = [mean(norm_gs, 'omitnan'), mean(norm_an, 'omitnan'), mean(norm_ps, 'omitnan')]; 
    
    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(mean_slopes);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
    %b.CData = [default_colors(7, :); default_colors(5, :); [0, 0, 0]; default_colors(6, :)];
    x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
    y = [norm_gs'; norm_an'; norm_ps'];
    plot(x, y, 'ko')
    hold off
    xticks([1, 2, 3])
    xticklabels(["Galvanic", "Anodic", "Pulsatile"])
    ylabel("Change in Firing Rate Slope (spk/s^2)")
    title("Recurrent Excitation Metric")
end