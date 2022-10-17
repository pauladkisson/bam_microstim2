%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Coefficient of Variation (CV)
function plot_cv(sim_name, sim_names, pulse_amps, stim_amps, t, N, top_N, num_group, ...
                 win_size, cv_window, default_colors, ex_brain, ex_c, ex_trial, ...
                 ex_neuron, brains, num_brains, pulse_coherences, galvanic_coherences, control_coherences, ...
                 start_trial, end_trial, num_trials, plot_name)
    if plot_name == "ex_neuron" || plot_name == "ex_trial"
        figure;
        hold on
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
            isi = zeros(length(t), N);
            for nn = 1:N
                spikeidx = recspikes(int2str(nn));
                if length(spikeidx) < 2
                    isi(:, nn) = NaN;
                    continue
                end
                neuron_isi_idx = diff(spikeidx);
                pre_isi_idx = max(neuron_isi_idx(1), spikeidx(1));
                isi(1:spikeidx(1)-1, nn) = t(pre_isi_idx);
                for i = 1:length(neuron_isi_idx)-1
                    isi_idx = neuron_isi_idx(i);
                    spikei = spikeidx(i);
                    next_spikei = spikeidx(i+1);
                    isi(spikei:next_spikei-1, nn) = t(isi_idx);
                end
                post_isi_idx = max(neuron_isi_idx(end), length(t)-spikeidx(end));
                isi(spikeidx(end):end, nn) = t(post_isi_idx);
            end
            p1_isi = isi(cv_window, 1:num_group);
            cv = movstd(p1_isi, win_size, [], 1) ./ movmean(p1_isi, win_size, 1);
            omit_idx = find(any(isnan(cv), 1));
            new_ex_neuron = ex_neuron;
            for idx = omit_idx
                if ex_neuron > omit_idx
                    new_ex_neuron = new_ex_neuron-1;
                elseif ex_neuron == omit_idx
                    disp("Example Neuron has less than 2 spikes, so it was omitted")
                    assert(1==2)
                end
            end
            cv(:, any(isnan(cv), 1)) = []; %omitnan
            if plot_name == "ex_neuron"
                plot(t(cv_window), cv(:, new_ex_neuron), 'Color', stim_color)
                scatter(t(recspikes(int2str(ex_neuron))), 10+j, 'Marker', '|', ...
                    "MarkerFaceColor", stim_color, "MarkerEdgeColor", stim_color)
            else
                plot(t(cv_window), mean(cv(:, 1:top_N), 2), 'Color', stim_color)
            end
        end
        hold off
        ylabel("P1 CV")
        xlabel("Time (s)")
        legend(["Pulsatile", "Galvanic", "Control"])
    
    elseif plot_name == "p1_wins"
        ball_rs = zeros(num_brains, num_group);
        figure(1);
        for sim_name = sim_names
            stim_cv = zeros(length(stim_amps), num_brains, num_group);
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
                    if ~contains(sim_name, "Discon")
                        num_wins = sum(final_decisions(:, stim_coherences==0)==1, 'all') * ones(num_group, 1);
                    else
                        num_wins = num_trials;
                        stim_coherences = 0;
                    end
                    c = 0;
                    for trial = start_trial:end_trial
                        relative_trial = trial - start_trial + 1;
                        if ~contains(sim_name, "Discon") && final_decisions(relative_trial, stim_coherences==c) ~= 1
                            continue %skip trials where P1 doesn't win
                        end
                        load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                            "recspikes")
                        isi = zeros(length(t), num_group);
                        for nn = 1:num_group
                            spikeidx = recspikes(int2str(nn));
                            if length(spikeidx) < 2
                                isi(:, nn) = NaN;
                                continue
                            end
                            neuron_isi_idx = diff(spikeidx);
                            pre_isi_idx = max(neuron_isi_idx(1), spikeidx(1));
                            isi(1:spikeidx(1)-1, nn) = t(pre_isi_idx);
                            for i = 1:length(neuron_isi_idx)-1
                                isi_idx = neuron_isi_idx(i);
                                spikei = spikeidx(i);
                                next_spikei = spikeidx(i+1);
                                isi(spikei:next_spikei-1, nn) = t(isi_idx);
                            end
                            post_isi_idx = max(neuron_isi_idx(end), length(t)-spikeidx(end));
                            isi(spikeidx(end):end, nn) = t(post_isi_idx);
                        end
                        cv = std(isi(cv_window, :), [], 1) ./ mean(isi(cv_window, :), 1);
                        nan_cv = isnan(cv);
                        cv = cv(~nan_cv);
                        num_wins(nan_cv) = num_wins(nan_cv) - 1;
                        stim_cv(j, brain, ~nan_cv) = reshape(stim_cv(j, brain, ~nan_cv), size(cv)) + cv;
                    end
                    if ~contains(sim_name, "Discon")
                        stim_cv(j, brain, :) = reshape(stim_cv(j, brain, :), size(num_wins)) ./ num_wins;
                    else
                        stim_cv(j, brain, :) = stim_cv(j, brain, :) / num_wins;
                    end
                end
            end
            pulse_cv = reshape(stim_cv(1, :, :), [num_brains, num_group]);
            galvanic_cv = reshape(stim_cv(2, :, :), [num_brains, num_group]);
            control_cv = reshape(stim_cv(3, 1, :), [1, num_group]);

            top_N = num_group * 0.2; %show top 20% closest to the electrode
            total_N = top_N*num_brains;
            tot_pulse_cv = reshape(pulse_cv(:, 1:top_N), [total_N, 1]);
            tot_galvanic_cv = reshape(galvanic_cv(:, 1:top_N), [total_N, 1]);
            tot_ctrl_cv = reshape(control_cv(1:top_N), [top_N, 1]);
            tot_ball_rs = reshape(ball_rs(:, 1:top_N), [total_N, 1]);
            tot_ctrl_ball_r = ctrl_ball_r(1:top_N);

            figure;
            hold on
            scatter(tot_ball_rs*1e6, tot_pulse_cv, [], default_colors(7, :), 'filled')
            scatter(tot_ball_rs*1e6, tot_galvanic_cv, [], default_colors(5, :), 'filled')
            scatter(tot_ctrl_ball_r*1e6, tot_ctrl_cv, [], "k", 'filled')
            hold off
            xlabel("Distance from Electrode (um)")
            ylabel("Coefficient of Variation (unitless)")
            title(strcat(sim_name, ": top 20% of P1"))

            top_N = num_group;
            total_N = top_N*num_brains;
            tot_pulse_cv = reshape(pulse_cv(:, 1:top_N), [total_N, 1]);
            tot_galvanic_cv = reshape(galvanic_cv(:, 1:top_N), [total_N, 1]);
            tot_ctrl_cv = reshape(control_cv(1:top_N), [top_N, 1]);
            tot_ball_rs = reshape(ball_rs(:, 1:top_N), [total_N, 1]);
            tot_ctrl_ball_r = ctrl_ball_r(1:top_N);

            figure;
            hold on
            scatter(tot_ball_rs*1e6, tot_pulse_cv, [], default_colors(7, :), 'filled')
            scatter(tot_ball_rs*1e6, tot_galvanic_cv, [], default_colors(5, :), 'filled')
            scatter(tot_ctrl_ball_r*1e6, tot_ctrl_cv, [], "k", 'filled')
            hold off
            xlabel("Distance from Electrode (um)")
            ylabel("Coefficient of Variation (unitless)")
            title(strcat(sim_name, ": Full P1"))

            if contains(sim_name, "Discon")
                figure(1);
                hold on
                scatter(tot_ball_rs*1e6, tot_pulse_cv, [], default_colors(7, :), 'filled')
                scatter(tot_ball_rs*1e6, tot_galvanic_cv, [], default_colors(5, :), 'filled')
                scatter(tot_ctrl_ball_r*1e6, tot_ctrl_cv, [], "k", 'filled')
                hold off
                xlabel("Distance from Electrode (um)")
                ylabel("Coefficient of Variation (unitless)")
            else
                figure(1);
                hold on
                scatter(tot_ball_rs*1e6, tot_pulse_cv, [], default_colors(7, :), 'filled', '^')
                scatter(tot_ball_rs*1e6, tot_galvanic_cv, [], default_colors(5, :), 'filled', '^')
                scatter(tot_ctrl_ball_r*1e6, tot_ctrl_cv, [], "k", 'filled', '^')
                hold off
                xlabel("Distance from Electrode (um)")
                ylabel("Coefficient of Variation (unitless)")
            end

            %trial/pop means
            pulse_trialpop_mean = mean(pulse_cv, 2);
            galvanic_trialpop_mean = mean(galvanic_cv, 2);
            ctrl_mean = mean(control_cv);
            pulse_norm_mean = pulse_trialpop_mean - ctrl_mean;
            galvanic_norm_mean = galvanic_trialpop_mean - ctrl_mean;
            stim_means = [mean(galvanic_norm_mean); mean(pulse_norm_mean)];
            stim_sems = [std(galvanic_norm_mean)/sqrt(num_brains); std(pulse_norm_mean)/sqrt(num_brains)];

            [~, p] = ttest2(galvanic_norm_mean, pulse_norm_mean)
            galvanic_norm_bar = mean(galvanic_norm_mean)
            galvanic_norm_sem = std(galvanic_norm_mean) ./ sqrt(length(galvanic_norm_mean))
            pulse_norm_bar = mean(pulse_norm_mean)
            pulse_norm_sem = std(pulse_norm_mean) ./ sqrt(length(pulse_norm_mean))

            if ~contains(sim_name, "Discon")
                connected_stim_means = stim_means;
                connected_gs_norm_mean = galvanic_norm_mean;
                connected_ps_norm_mean = pulse_norm_mean;
            else
                disconnected_stim_means = stim_means;
                disconnected_gs_norm_mean = galvanic_norm_mean;
                disconnected_ps_norm_mean = pulse_norm_mean;
            end
        end
        legend(["Connected", "", "", "Disconnected", "", ""])
        barcolor = [default_colors(5, :);
                    default_colors(7, :)];
        dis_pos = [1, 2];
        con_pos = [4, 5];
        figure;
        hold on
        b = bar(dis_pos, disconnected_stim_means, ...
            'FaceColor', 'flat', 'CData', barcolor);
        b = bar(con_pos, connected_stim_means, ...
            'FaceColor', 'flat', 'CData', barcolor);
        ylabel("Change in Coefficient of Variation (unitless)")
        plot([ones(1, num_brains)*dis_pos(1); ones(1, num_brains)*dis_pos(2)], [disconnected_gs_norm_mean'; disconnected_ps_norm_mean'], 'ko-')
        plot([ones(1, num_brains)*con_pos(1); ones(1, num_brains)*con_pos(2)], [connected_gs_norm_mean'; connected_ps_norm_mean'], 'ko-')
        hold off
        xticks([mean(dis_pos), mean(con_pos)])
        xticklabels(["Dis", "Con"])
        ylim([-0.0325, 0])

        [~, ~, stats] = anova1([connected_gs_norm_mean, disconnected_gs_norm_mean, connected_ps_norm_mean, disconnected_ps_norm_mean]);
        c = multcompare(stats)
    end
end