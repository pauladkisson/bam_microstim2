%%% Paul Adkisson
%%% 9/29/22
%%% Plot Oscillatory Modulataion of firing rates (mhat)
function plot_oscillation(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, ...
                        win_start, win_stop, default_colors, brains, num_brains, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        f0, start_trial, end_trial, num_trials, plot_name)
    T = win_stop - win_start;
    dt = t(2) - t(1);
    fs = 1/dt;
    nyquist = 1 / (2*dt); %Nyquist Frequency
    win_len = floor(T*fs);
    window = ones(win_len, 1);
    n_overlap = 0;
    stim_mhats = zeros(length(stim_amps), num_brains, num_group);
    ball_rs = zeros(num_brains, num_group);
    for brain = brains
        fprintf("brain %0.0f \n", brain)
        load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, brain]), "ball_r")
        ball_rs(brain, :) = ball_r;
        for j = 1:length(stim_amps)
            fprintf("j = %0.0f \n", j)
            stim_amp = stim_amps(j);
            pulse = j<= length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                    [sim_name, brain, stim_amp*1e9]);
                stim_coherences = pulse_coherences;
            else
                output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                    [sim_name, brain, stim_amp*1e9]);
                if stim_amp == 0
                    stim_coherences = control_coherences;
                    ctrl_ball_r = ball_r;
                    if brain ~= 1
                        continue
                    end
                else
                    stim_coherences = galvanic_coherences;
                end
            end
            load(strcat(output_stimpath, "/decisions.mat"), "final_decisions")
            if plot_name == "p1_wins"
                num_wins = sum(final_decisions==1, 'all');
                analyze_coherences = stim_coherences;
            elseif plot_name == "p1_loses"
                num_wins = sum(final_decisions==2, 'all');
                analyze_coherences = stim_coherences;
            elseif plot_name == "ex_c"
                num_wins = num_trials;
                analyze_coherences = ex_c;
            end
            for ex_c = analyze_coherences
                fprintf("Coherence %0.1f%% \n", ex_c*100)
                for trial = start_trial:end_trial
                    relative_trial = trial - start_trial + 1;
                    if (plot_name == "p1_wins" && final_decisions(relative_trial, stim_coherences==ex_c) ~= 1) || ...
                            (plot_name == "p1_loses" && final_decisions(relative_trial, stim_coherences==ex_c) ~= 2)
                        continue
                    end
                    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, trial])), ...
                        "recspikes")
                    spikes = zeros(length(t), num_group);
                    for nn = 1:num_group
                        for spike_idx = recspikes(int2str(nn))
                            spikes(spike_idx, nn) = 1;
                        end
                    end
                    %agg_spikes = sum(spikes(t>=win_start&t<win_stop, :), 2); %aggregate neurons for better resolution
                    for nn = 1:num_group
                        nn_spikes = spikes(t>=win_start&t<win_stop, nn);
                        [psd, freqs] = pwelch(nn_spikes, window, n_overlap, win_len);
                        freqs = freqs * nyquist / pi;
                        psd_snr = (psd - mean(psd(freqs>=100)) ) / std(psd(freqs>=100));
                        [~, f0_ind] = min(abs(freqs-f0));
                        peak_snr = psd_snr(f0_ind);
                        r0_hat = sum(nn_spikes) / T;
                        if peak_snr > 0 %peak SNR must be positive to give real m hat
                            mhat = get_mhat(peak_snr, r0_hat, T, fs, win_len);
                        else
                            mhat = 0;
                        end
                        stim_mhats(j, brain, nn) = stim_mhats(j, brain, nn) + mhat;
                    end
                end
            end
            stim_mhats(j, :, :) = stim_mhats(j, :, :) ./ num_wins;
        end
    end
    ctrl_mean = mean(stim_mhats(3, 1, :));
    ctrl_std = std(stim_mhats(3, 1, :));
    top_N = floor(num_group);
    total_N = top_N*num_brains;
    pulse_mhats = reshape(stim_mhats(1, :, 1:top_N), [total_N, 1]);
    galvanic_mhats = reshape(stim_mhats(2, :, 1:top_N), [total_N, 1]);
    ctrl_mhats = reshape(stim_mhats(3, 1, 1:top_N), [top_N, 1]);
    ball_rs = reshape(ball_rs(:, 1:top_N), [total_N, 1]);

    figure;
    set(gca, 'fontsize', 18)
    hold on
    scatter(ball_rs*1e6, pulse_mhats, [], ones(total_N, 3).*default_colors(7, :), 'filled')
    scatter(ball_rs*1e6, galvanic_mhats, [], ones(total_N, 3).*default_colors(5, :), 'filled')
    scatter(ball_rs*1e6, ctrl_mhats, [], "k", 'filled')
    yline(ctrl_mean, 'k--')
    hold off
    xlabel("Distance from Electrode ($\mu$m)", "Interpreter", "LaTeX")
    ylabel("Modulation Strength $\hat{m}$", 'Interpreter', 'LaTeX')
    if sim_name == "EMBC Disconnected" || sim_name == "DepolBlockDiscon" || sim_name == "Test"
        title("Disconnected")
    else
        title("Connected")
    end
    
    %trial/pop means
    pulse_trialpop_mean = mean(reshape(stim_mhats(1, :, :), [num_brains, num_group]), 2);
    galvanic_trialpop_mean = mean(reshape(stim_mhats(2, :, :), [num_brains, num_group]), 2);
    pulse_norm_mean = pulse_trialpop_mean - ctrl_mean;
    galvanic_norm_mean = galvanic_trialpop_mean - ctrl_mean;
    stim_means = [mean(galvanic_norm_mean); mean(pulse_norm_mean)];
    stim_sems = [std(galvanic_norm_mean)/sqrt(num_brains); std(pulse_norm_mean)/sqrt(num_brains)];

    figure;
    set(gca, 'fontsize', 18)
    hold on
    b = bar(stim_means);
    b.FaceColor = 'flat';
    b.CData = [default_colors(5, :); default_colors(7, :)];
    plot([ones(1, num_brains); 2*ones(1, num_brains)], [galvanic_norm_mean'; pulse_norm_mean'], 'ko-')
    hold off
    xticks([1, 2])
    xticklabels(["Galvanic", "Pulsatile"])
    ylabel("Change in Modulation Index $\hat{m}$", 'Interpreter', 'LaTeX')
    if sim_name == "EMBC Disconnected" || sim_name == "DepolBlockDiscon" || sim_name == "Test"
        title("Disconnected")
    else
        title("Connected")
    end
    
end


function mhat = get_mhat(peak_snr, r0, T, fs, win_len)
    N_wins = ceil(T*fs / win_len);
    peak_snr_cor = peak_snr*sqrt(N_wins); %welch correction
    mhat = 2*sqrt(peak_snr_cor / (r0*T));
end