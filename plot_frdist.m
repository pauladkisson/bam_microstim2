%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Firing Rate Distribution over distance from electrode
function plot_frdist(sim_names, ex_c, pulse_amps, stim_amps, t, t_cut, num_group, num_affected, ...
                     win_start, win_stop, default_colors, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     anodic_coherences, start_trial, end_trial, num_trials, plot_name)
    for sim_name = sim_names
        stim_frs = zeros(length(stim_amps), num_trials, num_group);
        load(sprintf("Simulation %s/ustim/r.mat", sim_name), "ball_r")
        ball_rs = get_ball_rs(ball_r, num_affected, num_group);
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            c = ex_c(j);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
                stim_coherences = pulse_coherences;
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
                if stim_amp < 0 %cathodic GS
                   stim_coherences = galvanic_coherences;
                elseif stim_amp == 0
                    stim_coherences = control_coherences;
                else %anodic GS
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

        %Full Population Aggregated Activity
        popmean_pulse = reshape(mean(stim_frs(1, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_galvanic = reshape(mean(stim_frs(2, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_ctrl = reshape(mean(stim_frs(3, :, :), 3, 'omitnan'), [num_trials, 1]);
        popmean_anodic = reshape(mean(stim_frs(4, :, :), 3, 'omitnan'), [num_trials, 1]);
        norm_pulse = popmean_pulse - popmean_ctrl;
        norm_galvanic = popmean_galvanic - popmean_ctrl;
        norm_anodic = popmean_anodic - popmean_ctrl;
        stim_means = [mean(norm_galvanic, 'omitnan'), mean(norm_anodic, 'omitnan'), ...
                      mean(norm_pulse, 'omitnan')];
        stim_stds = [std(norm_galvanic, [], 'omitnan'), std(norm_anodic, [], 'omitnan'), ...
                     std(norm_pulse, [], 'omitnan')];
        figure;
        set(gca, 'fontsize', 18)
        hold on
        b = bar(stim_means, 1);
        b.FaceColor = 'flat';
        b.CData = [default_colors(5, :); default_colors(6, :); default_colors(7, :)];
        x = [1, 2, 3];
        errorbar(x, stim_means, stim_stds, 'k.', 'Linewidth', 20, 'Capsize', 0)
        hold off
        xticks([1, 2, 3])
        xticklabels(["Galvanic", "Anodic", "Pulsatile"])
        ylabel("Change in Firing Rate (spk/s)")
        ylim([-4, 4])
        title(sim_name)
        
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
end