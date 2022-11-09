%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current for a variety of PRs and
%%% Amplitudes
function GeneratePopMicroStim(t, t_task, t_taskoff, stim_duration, stim_freqs, ...
            gL, stim_amps, save_amp, N, sim_path, plot_ustim)
    I_ustim_base = zeros(length(t), N);
    dt = t(2) - t(1);
    if ~isempty(stim_freqs)
        for j = 1:length(stim_freqs)
            stim_amp = stim_amps(j);
            stim_freq = stim_freqs(j);
            for i = 1:length(t) 
                if floor(1/(stim_freq*dt)) >= length(t)
                    mod_i = i;
                else
                    mod_i = mod(i, floor(1/(stim_freq*dt))) + 1;
                end
                if t(mod_i) <=  stim_duration
                    I_ustim_base(i, j) = stim_amp;
                elseif t(mod_i) > stim_duration && t(mod_i) <= 2*stim_duration
                    I_ustim_base(i, j) = -stim_amp;
                end
            end
        end
    else
        I_ustim_base = zeros(length(t), N);
        for j = 1:length(stim_amps)
            stim_amp = stim_amps(j);
            I_ustim_base(:, j) = stim_amp;
        end
    end

    I_ustim_base(t<t_task|t>t_taskoff, :) = 0;

    z_thia = 0.002; %Thia's z=2mm
    ball_r = ones(1, N)*z_thia;
    electric_r = mirror_est(ball_r);
    I_ustim = I_ustim_base.*electric_r;
    true_amps = I_ustim(t==t_task, :);
    Vmir = true_amps ./ gL;

    basepath = strcat(sim_path, "/ustim");
    mkdir(basepath)
    save(strcat(basepath, "/r.mat"), 'ball_r')
    if ~isempty(stim_freqs)
        save(strcat(basepath, sprintf("/%0.2fuA_pulse.mat", save_amp*1e6)), "I_ustim", 'Vmir')
    else
        save(strcat(basepath, sprintf("/%0.2fuA_galvanic.mat", save_amp*1e6)), "I_ustim", 'Vmir')
    end

    if plot_ustim
        figure;
        hold on
        scatter(ball_r*1e6, Vmir*1e3)
        xlabel("Distance from Electrode (um)")
        ylabel("Stimulation Depolarization (mV)")

        figure;
        scatter(ball_r*1e6, Vmir*gL*1e9)
        xlabel("Distance from Electrode (um)")
        ylabel("Internal Stimulation Amplitude (nA)")
    end
end