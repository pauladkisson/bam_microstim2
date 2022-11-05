%%% Paul Adkisson
%%% 2/14/2022
%%% Plot example Population Firing Rates
function plot_frs(sim_name, pulse_amps, stim_amps, p, f, N, N_E, t, t_task,...
                  t_taskoff, default_colors, ex_stim_j, ex_c, ex_trial, plot_name)
    dt = t(2) - t(1);
    if plot_name == "single_stim"
        pulse = ex_stim_j<=length(pulse_amps);
        stim_amp = stim_amps(ex_stim_j);
        if pulse
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                [sim_name, stim_amp*1e6]);
        elseif stim_amp == 0
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        else
            output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                [sim_name, stim_amp*1e6]);
        end
        load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "recspikes")
        [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
        figure;
        set(gca, 'fontsize', 18)
        hold on
        for i = 1:p+2
            plot(t, pop_frs(:, i), 'Linewidth', 4)
        end
        xline(t_task, 'k--')
        xline(t_taskoff, 'k--')
        hold off
        xlabel("Time (s)")
        ylabel("Population Firing Rate (spk/s)")
        legend(["P1", "P2", "NS", "Int"], ...
            "Location", "northwest")


    elseif plot_name == "subplot"
        figure;
        set(gca, 'fontsize', 18)
        axs = zeros(length(stim_amps), 1);
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
            elseif stim_amp == 0
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "recspikes")
            [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
            mean_fr = mean(pop_frs(t>=2.5&t<3, 1));
            axs(j) = subplot(length(stim_amps), 1, j);
            set(axs(j), 'fontsize', 18);
            hold on
            plot(t, pop_frs, 'Linewidth', 4)
            hold off
            if pulse
                title("Pulsatile Stimulation")
                fprintf("Pulse mean fr: %0.3f \n", mean_fr) 
            elseif stim_amp == 0
                title("Control")
                fprintf("Control mean fr: %0.3f \n", mean_fr) 
            else
                title("Galvanic Stimualtion")
                fprintf("Galvanic mean fr: %0.3f \n", mean_fr) 
            end
            if j == length(stim_amps)
                xlabel("Time (s)")
            elseif j == ceil(length(stim_amps)/2)
                ylabel("Population Firing Rate (spk/s)")
            end
        end
        linkaxes(axs)

    elseif plot_name == "p1_only"
        figure;
        set(gca, 'fontsize', 18)
        hold on
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            pulse = j<=length(pulse_amps);
            if pulse
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_pulse", ...
                    [sim_name, stim_amp*1e6]);
            elseif stim_amp == 0
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
            else
                output_stimpath = sprintf("Simulation %s/data/%0.2fuA_galvanic", ...
                    [sim_name, stim_amp*1e6]);
            end
            load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "recspikes")
            [pop_frs, ~] = recspikes2popfrs(recspikes, t, N, dt, p, f, N_E);
            if pulse
                plot(t, pop_frs(:, 1), 'Color', default_colors(7, :), 'Linewidth', 2)
            elseif stim_amp < 0 %cathodic gs
                plot(t, pop_frs(:, 1), 'Color', default_colors(5, :), 'Linewidth', 2)
            elseif stim_amp == 0
                plot(t, pop_frs(:, 1), "k", 'Linewidth', 2)
            else %anodic gs
                plot(t, pop_frs(:, 1), 'Color', default_colors(6, :), 'Linewidth', 2)
            end
        end
        xlabel("Time (s)")
        ylabel("P1 Firing Rate (spk/s)")
    end
end