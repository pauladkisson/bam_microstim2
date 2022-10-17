%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current
function GenerateMicroStim_Desync(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                           pulse_amps, dc_amps, N, num_group, brains, ...
                           sim_path)
    I_ustim_base = zeros(length(t), N);
    dt = t(2) - t(1);
    stim_amps = [pulse_amps, dc_amps];
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        is_pulse = j <= length(pulse_amps);
        if is_pulse
            desync_jitter = randi(floor(1/(stim_freq*dt)), N, 1);
            for i = 1:length(t) 
                if floor(1/(stim_freq*dt)) >= length(t)
                    mod_i = i*ones(N, 1);
                else
                    mod_i = mod(i+desync_jitter, floor(1/(stim_freq*dt))) + 1;
                end
                if any(t(mod_i) <=  stim_duration)
                    I_ustim_base(i, t(mod_i)<=stim_duration) = stim_amp;
                elseif any(t(mod_i) > stim_duration & t(mod_i) <= 2*stim_duration)
                    I_ustim_base(i, t(mod_i) > stim_duration & t(mod_i) <= 2*stim_duration) = -stim_amp;
                end
            end
        else
            I_ustim_base = ones(length(t), N)*stim_amp;
        end
        I_ustim_base(t<t_task|t>t_taskoff, :) = 0;
        for brain = brains
            brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
            load(strcat(brainpath, "/r.mat"), "electric_r")
            I_ustim = [I_ustim_base(:, 1:num_group).*electric_r, zeros(length(t), N-num_group)];
            
            if is_pulse
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_pulse.mat", stim_amp*1e9)), "I_ustim")
            else
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_galvanic.mat", stim_amp*1e9)), "I_ustim")
            end
        end
        %{
        figure;
        plot(t, I_ustim(:, 1:num_group)*1e9)
        title(sprintf("Stim Amp: %0.0fnA", stim_amp*1e9))
        
        figure;
        scatter(1:N, jitter_mod_i)
        %}
    end
end