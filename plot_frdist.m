%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_names, ex_c, pulse_amps, stim_amps, t, t_cut, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     anodic_coherences, start_trial, end_trial, num_trials, plot_name)
    sim_pulse = zeros(length(sim_names), num_trials);
    sim_galvanic = zeros(length(sim_names), num_trials);
    sim_anodic = zeros(length(sim_names), num_trials);
    for sim_name = sim_names
        disp(" ")
        disp(sim_name)
        stim_frs = zeros(length(stim_amps), num_trials, num_group);
        load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
        ball_rs = get_ball_rs(ball_r, num_affected, num_group);
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            c = ex_c(j);
            if pulse
                disp("Pulse")
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
                stim_coherences = pulse_coherences;
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
                if stim_amp < 0 %cathodic GS
                    disp("Cathodic GS")
                   stim_coherences = galvanic_coherences;
                elseif stim_amp == 0
                    disp("Control")
                    stim_coherences = control_coherences;
                else %anodic GS
                    disp("Anodic GS")
                    stim_coherences = anodic_coherences;
                end
            end
            try
                load(strcat(output_stimpath, "/decisions.mat"), "decisions", "decision_times")
                control_decs = load(sprintf("Simulation %s/data/0.00uA_galvanic/decisions.mat", sim_name));
                ctrl_decs = control_decs.decisions;
            catch
                assert(sim_name~="Brainless_m=0_Con")
            end
            for trial = start_trial:end_trial
                relative_trial = trial - start_trial + 1;
                if (plot_name == "p1_wins" && (decisions(relative_trial, stim_coherences==c) ~= 1 || ...
                        ctrl_decs(relative_trial, control_coherences==ex_c(3))~=1)) || ...
                   (plot_name == "p1_loses" && ( decisions(relative_trial, stim_coherences==c) ~= 2 || ...
                        ctrl_decs(relative_trial, control_coherences==ex_c(3))~=2)) || ...
                   (contains(plot_name, "p1") && decision_times(relative_trial, stim_coherences==c) > t_cut)
                    %skip trials P1 doesn't win/lose for stim and control
                    %or decision takes too long
                    stim_frs(j, relative_trial, :) = NaN;
                    continue
                end
                load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                    "recspikes")
                g1_taskfrs = zeros(num_group, 1);
                for nn = 1:num_group
                    spiketimes = t(recspikes(int2str(nn)));
                    g1_taskfrs(nn) = sum(spiketimes>=win_start & spiketimes<win_stop) / (win_stop - win_start);
                end
                stim_frs(j, trial, :) = g1_taskfrs;
            end
        end
        pulse_frs = reshape(mean(stim_frs(1, :, :), 2, 'omitnan'), [num_group, 1]);
        galvanic_frs = reshape(mean(stim_frs(2, :, :), 2, 'omitnan'), [num_group, 1]);
        control_frs = reshape(mean(stim_frs(3, :, :), 2, 'omitnan'), [num_group, 1]);
        anodic_frs = reshape(mean(stim_frs(4, :, :), 2, 'omitnan'), [num_group, 1]);

        pulse_sems = reshape(std(stim_frs(1, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
        galvanic_sems = reshape(std(stim_frs(2, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
        control_sems = reshape(std(stim_frs(3, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);
        anodic_sems = reshape(std(stim_frs(4, :, :), [], 2, 'omitnan'), [num_group, 1]) / sqrt(num_trials);

        figure;
        set(gca, 'fontsize', 18)
        hold on
        errorbar(ball_rs*1e6, pulse_frs, pulse_sems, '.', ...
            'Color', default_colors(7, :), 'MarkerSize', 20)
        errorbar(ball_rs*1e6, galvanic_frs, galvanic_sems, '.', ...
            'Color', default_colors(5, :), 'MarkerSize', 20)
        errorbar(ball_rs*1e6, control_frs, control_sems, "k.", 'MarkerSize', 20)
        errorbar(ball_rs*1e6, anodic_frs, anodic_sems, '.', ...,
            'MarkerSize', 20, 'Color', default_colors(6, :))
        if plot_name == "ex_c_<400"
            xlim([0, 400])
        elseif plot_name == "ex_c_>400"
            xlim([400, 2000])
        elseif plot_name == "ex_c_zoom"
            if contains(sim_name, "Discon")
                ylim([15, 40])
            elseif contains(sim_name, "Int")
                ylim([0, 1])
            elseif contains(sim_name, "Rec")
                ylim([40, 70])
            end
            xlim([0, 2000])
        end
        hold off
        xlabel("Distance from Electrode (um)")
        ylabel("Firing Rate (spk/s)")
        title(sim_name)
        
        %Statistics
        pulse_frs = pulse_frs(1:num_affected);
        galvanic_frs = galvanic_frs(1:num_affected);
        ball_rs = ball_rs(1:num_affected);
        [~, p_ps_cgs] = kstest2(pulse_frs, galvanic_frs);
        p_ps = ones(num_affected, 1);
        p_gs = ones(num_affected, 1);
        p_an = ones(num_affected, 1);
        p_cgs_ags = ones(num_affected, 1);
        p_ps_less = ones(num_affected, 1);
        p_gs_less = ones(num_affected, 1);
        for neuron = 1:num_affected
            ps_frs = reshape(stim_frs(1, :, neuron), [num_trials, 1]);
            gs_frs = reshape(stim_frs(2, :, neuron), [num_trials, 1]);
            ctrl_frs = reshape(stim_frs(3, :, neuron), [num_trials, 1]);
            an_frs = reshape(stim_frs(4, :, neuron), [num_trials, 1]);
            [~, p_ps(neuron)] = ttest2(ps_frs, ctrl_frs, 'Tail', 'right');
            [~, p_gs(neuron)] = ttest2(gs_frs, ctrl_frs, 'Tail', 'right');
            [~, p_an(neuron)] = ttest2(an_frs, ctrl_frs, 'Tail', 'left');
            [~, p_cgs_ags(neuron)] = ttest2(abs(gs_frs-ctrl_frs), abs(an_frs-ctrl_frs));
            [~, p_ps_less(neuron)] = ttest2(ps_frs, ctrl_frs, 'Tail', 'left');
            [~, p_gs_less(neuron)] = ttest2(gs_frs, ctrl_frs, 'Tail', 'left');
        end
        % Bonferroni correction
        p_ps = p_ps * num_affected;
        p_gs = p_gs * num_affected;
        p_an = p_an * num_affected;
        p_cgs_ags = p_cgs_ags * num_affected;
        p_ps_less = p_ps_less * num_affected;
        p_gs_less = p_gs_less * num_affected;
        sig_thresh = 0.05;
        ps_sig_dist = max(ball_rs(p_ps < sig_thresh));
        gs_sig_dist = max(ball_rs(p_gs < sig_thresh));
        an_sig_dist = max(ball_rs(p_an < sig_thresh));
        cgs_ags_dist = max(ball_rs(p_cgs_ags < sig_thresh));
        ps_num_less = sum(p_ps_less < sig_thresh);
        gs_num_less = sum(p_gs_less < sig_thresh);
        disp("DISTRIBUTIONS") 
        fprintf([...
            'PS and CGS induce different firing rate distributions ', ...
            '(p=%0.1e). \n'], p_ps_cgs)
        fprintf('CGS induces significant FR increases up to %0.1fum away. \n', ...
                gs_sig_dist*1e6)
        fprintf('PS induces significant FR increases up tp %0.1fum away. \n', ...
                ps_sig_dist*1e6)
        fprintf('AGS induces significant FR decreases up tp %0.1fum away. \n', ...
                an_sig_dist*1e6)
        fprintf([...
            'CGS and AGS had equal and opposite effects for neurons >%0.1fum ', ...
            'away. \n'], cgs_ags_dist*1e6)
        fprintf('PS had %0.0f neurons less than control. \n', ps_num_less)
        fprintf('GS had %0.0f neurons less than control. \n', gs_num_less)

        %Full Population Aggregated Activity
        popmean_pulse = reshape(mean(stim_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_galvanic = reshape(mean(stim_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_ctrl = reshape(mean(stim_frs(3, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_anodic = reshape(mean(stim_frs(4, :, :), 3, 'omitnan'), [num_trials, 1]);
        mean_ctrl = mean(popmean_ctrl);
        norm_pulse = popmean_pulse - mean_ctrl;
        norm_galvanic = popmean_galvanic - mean_ctrl;
        norm_anodic = popmean_anodic - mean_ctrl;
        norm_control = popmean_ctrl - mean_ctrl;
        ps_quantiles = quantile(norm_pulse, [0.25, 0.5, 0.75]);
        cgs_quantiles = quantile(norm_galvanic, [0.25, 0.5, 0.75]);
        ags_quantiles = quantile(norm_anodic, [0.25, 0.5, 0.75]);
        ctrl_quantiles = quantile(norm_control, [0.25, 0.5, 0.75]);
        
%         stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
%                       mean(norm_pulse, 'omitnan')];
%         stim_stds = [std(norm_galvanic, [], 'omitnan'), std(norm_anodic, [], 'omitnan'), ...
%                      std(norm_pulse, [], 'omitnan')];
%         stim_trials = [sum(~isnan(norm_galvanic)), sum(~isnan(norm_anodic)), ...
%                        sum(~isnan(norm_pulse))];
%         figure;
%         set(gca, 'fontsize', 18)
%         hold on
%         b = bar(stim_means, 1);
%         b.FaceColor = 'flat';
%         b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
%         x = [1, 2, 3];
%         errorbar(x, stim_means, stim_stds, 'k.', 'Linewidth', 20, 'Capsize', 0)
%         hold off
%         xticks([1, 2, 3])
%         xticklabels(["Galvanic", "Anodic", "Pulsatile"])
%         ylabel("Change in Firing Rate (spk/s)")
%         %ylim([-4, 4])
%         title(sim_name)
        figure;
        set(gca, 'fontsize', 18)
        hold on
        colors = [default_colors(5, :); default_colors(6, :); default_colors(7, :); [0, 0, 0]];
        boxplot([norm_galvanic, norm_anodic, norm_pulse, norm_control], 'PlotStyle', 'traditional', ...
            'Colors', colors, 'Symbol', ".")
        hold off
        xticks([1, 2, 3, 4])
        xticklabels(["Galvanic", "Anodic", "Pulsatile", "Control"])
        ylabel("Change in Firing Rate (spk/s)")
        ylim([-6, 6])
        title(sim_name)
        
        %Statistics
        disp("POPULATION AVERAGES")
        sim_pulse(sim_names==sim_name, :) = norm_pulse;
        sim_galvanic(sim_names==sim_name, :) = norm_galvanic;
        sim_anodic(sim_names==sim_name, :) = norm_anodic;
%         [~, p_cgs_ags] = ttest2(abs(norm_galvanic), abs(norm_anodic));
%         [~, p_ps_cgs] = ttest2(norm_pulse, norm_galvanic);
%         fprintf([...
%             'AGS induced a different average change in FR (%0.2f +/- %0.2fspk/s)', ...
%             ' than CGS (%0.2f +/- %0.3fspk/s), p=%0.1e. \n'], stim_means(2), ...
%             stim_stds(2)/sqrt(stim_trials(2)), stim_means(1), ...
%             stim_stds(1)/sqrt(stim_trials(1)), p_cgs_ags)
%         fprintf([...
%             'PS induced a different average change in FR (%0.2f +/- %0.2fspk/s)', ...
%             ' than CGS (%0.2f +/- %0.3fspk/s), p=%0.1e. \n'], stim_means(3), ...
%             stim_stds(3)/sqrt(stim_trials(3)), stim_means(1), ...
%             stim_stds(1)/sqrt(stim_trials(1)), p_ps_cgs)
        [p_median, ~, stats] = kruskalwallis([norm_pulse, norm_galvanic, ...
                                              norm_control, norm_anodic]);
        fprintf([...
            "Stimulation induces significantly different firing rates (p=%0.1e). \n"], ...
            p_median)
        c = multcompare(stats);
        p_ps_cgs = c(1, end);
        p_cgs_ags = ranksum(abs(norm_galvanic), abs(norm_anodic));
        fprintf([...
            'AGS induced a different average change in FR (%0.2f, %0.2f, %0.2f)', ...
            ' than CGS (%0.2f, %0.2f, %0.2f), p=%0.1e. \n'], ...
            ags_quantiles(1), ags_quantiles(2), ags_quantiles(3), ...
            cgs_quantiles(1), cgs_quantiles(2), cgs_quantiles(3), p_cgs_ags)
        fprintf([...
            'PS induced a different average change in FR (%0.2f +/- %0.2fspk/s)', ...
            ' than CGS (%0.2f +/- %0.3fspk/s), p=%0.1e. \n'],  ...
            ps_quantiles(1), ps_quantiles(2), ps_quantiles(3), ...
            cgs_quantiles(1), cgs_quantiles(2), cgs_quantiles(3), p_ps_cgs)
        
        %Unaffected P1
        popmean_pulse = reshape(mean(stim_frs(1, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
        popmean_galvanic = reshape(mean(stim_frs(2, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
        popmean_ctrl = reshape(mean(stim_frs(3, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
        popmean_anodic = reshape(mean(stim_frs(4, :, num_affected+1:end), 3, 'omitnan'), [num_trials, 1]);
        norm_pulse = popmean_pulse - popmean_ctrl;
        norm_galvanic = popmean_galvanic - popmean_ctrl;
        norm_anodic = popmean_anodic - popmean_ctrl;
        stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
                      mean(norm_pulse, 'omitnan')]; 
        figure;
        set(gca, 'fontsize', 18)
        hold on
        b = bar(stim_means);
        b.FaceColor = 'flat';
        b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
        x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
        y = [norm_galvanic'; norm_anodic'; norm_pulse'];
        plot(x, y, 'ko')
        hold off
        xticks([1, 2, 3])
        xticklabels(["Galvanic", "Anodic", "Pulsatile"])
        ylabel("Change in Firing Rate (spk/s)")
        %ylim([-4, 4])
        title("P1 Unaffected")
        
        %Affected P1 only
        popmean_pulse = reshape(mean(stim_frs(1, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
        popmean_galvanic = reshape(mean(stim_frs(2, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
        popmean_ctrl = reshape(mean(stim_frs(3, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
        popmean_anodic = reshape(mean(stim_frs(4, :, 1:num_affected), 3, 'omitnan'), [num_trials, 1]);
        norm_pulse = popmean_pulse - popmean_ctrl;
        norm_galvanic = popmean_galvanic - popmean_ctrl;
        norm_anodic = popmean_anodic - popmean_ctrl;
        stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
                      mean(norm_pulse, 'omitnan')];
        figure;
        set(gca, 'fontsize', 18)
        hold on
        b = bar(stim_means);
        b.FaceColor = 'flat';
        b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
        x = [ones(1, num_trials); 2*ones(1, num_trials); 3*ones(1, num_trials)];
        y = [norm_galvanic'; norm_anodic'; norm_pulse'];
        plot(x, y, 'ko')
        hold off
        xticks([1, 2, 3])
        xticklabels(["Galvanic", "Anodic", "Pulsatile"])
        ylabel("Change in Firing Rate (spk/s)")
        %ylim([-4, 4])
        title("Affected P1")
    end
    % Simulation-wise comparisons
    if length(sim_names) == 3 % discon, p1_int, p1_rec 
        disp(" ")
        disp("SIMULATION-WISE")
        p_ps_discon_inh = ranksum(sim_pulse(1, :), sim_pulse(2, :));
        p_gs_discon_inh = ranksum(sim_galvanic(1, :), sim_galvanic(2, :));
        p_ps_discon_rec = ranksum(sim_pulse(1, :), sim_pulse(3, :));
        p_gs_discon_rec = ranksum(sim_galvanic(1, :), sim_galvanic(3, :));
        p_an_discon_rec = ranksum(sim_anodic(1, :), sim_anodic(3, :));
        fprintf([...
                'PS (p=%0.1e) and CGS (%0.1e) were less effective under ', ...
                'feedback inhibition than disconnected. \n'], ...
                p_ps_discon_inh, p_gs_discon_inh)
        fprintf([...
                'PS (p=%0.1e) and CGS (%0.1e) were less/more effective under ', ...
                'recurrent excitation than disconnected respectively. \n'], ...
                p_ps_discon_rec, p_gs_discon_rec)
        fprintf([...
                'AGS (%0.1e) was more effective under ', ...
                'recurrent excitation than disconnected. \n'], ...
                p_an_discon_rec)
    end
end