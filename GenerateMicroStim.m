%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current
function GenerateMicroStim(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                          min_r, max_r, num_affected, thresh_cor, gL, ...
                          pulse_amps, dc_amps, N, sim_path, plot)
    I_ustim_base = zeros(length(t), N);
    dt = t(2) - t(1);
    stim_amps = [pulse_amps, dc_amps];
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        is_pulse = j <= length(pulse_amps);
        if is_pulse
            for i = 1:length(t) 
                if floor(1/(stim_freq*dt)) >= length(t)
                    mod_i = i;
                else
                    mod_i = mod(i, floor(1/(stim_freq*dt))) + 1;
                end
                if t(mod_i) <=  stim_duration
                    I_ustim_base(i, :) = stim_amp;
                elseif t(mod_i) > stim_duration && t(mod_i) <= 2*stim_duration
                    I_ustim_base(i, :) = -stim_amp;
                end
            end
        else
            I_ustim_base = ones(length(t), N)*stim_amp;
        end
        I_ustim_base(t<t_task|t>t_taskoff, :) = 0;
        
        regular_r = (max_r - min_r) / (num_affected-1);
        ball_r = min_r:regular_r:max_r;
        electric_r = mirror_est(ball_r);
        I_ustim = [I_ustim_base(:, 1:num_affected).*electric_r, zeros(length(t), N-num_affected)];
        true_amps = I_ustim(t==t_task, 1:num_affected);
        Vmir = true_amps ./ gL;

        basepath = strcat(sim_path, "/ustim");
        mkdir(basepath)
        save(strcat(basepath, "/r.mat"), 'ball_r')
        if is_pulse
            I_ustim = I_ustim * thresh_cor;
            Vmir = Vmir * thresh_cor;
            save(strcat(basepath, sprintf("/%0.2fuA_pulse.mat", stim_amp*1e6)), "I_ustim", 'Vmir')
        else
            save(strcat(basepath, sprintf("/%0.2fuA_galvanic.mat", stim_amp*1e6)), "I_ustim", 'Vmir')
        end
        if plot
            figure;
            hold on
            scatter(ball_r*1e6, Vmir*1e3)
            xlabel("Distance from Electrode (um)")
            ylabel("Stimulation Depolarization (mV)")
            if is_pulse
                title("Pulse")
            elseif stim_amp == 0
                title("Control")
            else
                title("Galvanic")
            end

            figure;
            scatter(ball_r*1e6, Vmir*gL*1e9)
            xlabel("Distance from Electrode (um)")
            ylabel("Internal Stimulation Amplitude (nA)")
            if is_pulse
                title("Pulse")
            elseif stim_amp == 0
                title("Control")
            else
                title("Galvanic")
            end
        end
    end
end