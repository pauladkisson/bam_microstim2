%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current
function GenerateMicroStim(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                          pulse_amps, dc_amps, N, num_group, brains, sim_path)
    thresh_cor = 0.211; %Threshold correction factor for pulses
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
        for brain = brains
            brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
            load(strcat(brainpath, "/r.mat"), "electric_r", "ball_r")
            I_ustim = [I_ustim_base(:, 1:num_group).*electric_r, zeros(length(t), N-num_group)];
            I_ustim(:, 2:2:end) = 0; %Every other neuron is unaffected
            gL = 25*1e-9;
            true_amps = I_ustim(t==t_task, 1:num_group);
            Vmir = true_amps ./ gL;
            
            if is_pulse
                I_ustim = I_ustim * thresh_cor;
                Vmir = Vmir * thresh_cor;
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_pulse.mat", stim_amp*1e9)), "I_ustim", 'Vmir')
            else
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_galvanic.mat", stim_amp*1e9)), "I_ustim", 'Vmir')
            end
        end
        %{
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
        %}
    end
end