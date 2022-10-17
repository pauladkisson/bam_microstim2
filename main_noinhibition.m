%%% Paul Adkisson
%%% 11.7.21
%%% Main Body of Wang (2002) Biophysical Attractor Model
%%% Adds Pulse-refractory period

clear;
sim_name = "NoInhibitionNoTaskTrials";
sim_path = sprintf("Simulation %s", sim_name);
tic;
load(strcat(sim_path, "/bam_constants.mat"))
load(strcat(sim_path, "/adja.mat"))
load(strcat(sim_path, "/conductances.mat"))
for brain = brains
    fprintf("Brain %0.0f \n", brain)
    brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            fprintf("Pulse Stimulation Amplitude: %0.1fnA \n", stim_amp*1e9)
            input_stimpath = strcat(brainpath, sprintf("/ustim/%0.1fnA_pulse.mat", stim_amp*1e9));
            output_stimpath = strcat(brainpath, sprintf("/data/%0.1fnA_pulse", stim_amp*1e9));
            stim_coherences = pulse_coherences;
        elseif stim_amp ~= 0
            fprintf("Galvanic Stimulation Amplitude: %0.1fnA \n", stim_amp*1e9)
            input_stimpath = strcat(brainpath, sprintf("/ustim/%0.1fnA_galvanic.mat", stim_amp*1e9));
            output_stimpath = strcat(brainpath, sprintf("/data/%0.1fnA_galvanic", stim_amp*1e9));
            stim_coherences = galvanic_coherences;
        else
            disp("0nA Stimulation Control")
            input_stimpath = strcat(brainpath, sprintf("/ustim/%0.1fnA_galvanic.mat", stim_amp*1e9));
            output_stimpath = strcat(brainpath, sprintf("/data/%0.1fnA_galvanic", stim_amp*1e9));
            stim_coherences = control_coherences;
        end
        load(input_stimpath)
        mkdir(output_stimpath)
        for c = stim_coherences
            fprintf("Coherence: %0.1f%% \n", c*100)
            input_coherentpath = sprintf("Simulation %s/spikes/c=%0.3f", sim_name, c);
            output_coherentpath = strcat(output_stimpath, sprintf("/c=%0.3f", c));
            mkdir(output_coherentpath)
            %parfor trial = start_trial:end_trial
            for trial = start_trial:end_trial
                fprintf("trial: %0.0f \n", trial)
                input_trialpath = strcat(input_coherentpath, sprintf("/trial%0.0f/input.mat", trial));
                output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));

                %external input
                input_spikes = load(input_trialpath, 'spikes');
                spikes = input_spikes.spikes;
                
                %ustim
                I_ustim_temp = I_ustim;

                %LIF variables
                RP_ind = zeros(1, N);
                if pulse
                    RP_pp_ind = zeros(1, N);
                    RP_ps_ind = zeros(1, N);
                end
                Vm = EL*ones(length(t), N);
                s_ampa = zeros(length(t), N);
                s_ampa_ext = zeros(length(t), N);
                s_nmda = zeros(length(t), N);
                x_nmda = zeros(length(t), N);
                s_gaba = zeros(length(t), N);

                for i = 1:length(t)-1
                    % Synaptic current
                    if i <= delay_ind
                        I_ch = zeros(1, N);
                    else
                        V_ch = Vm(i, :);
                        s_ampa_ch = s_ampa(i-delay_ind, :);
                        s_ampa_ext_ch = s_ampa_ext(i, :);
                        g_ampa_ext = s_ampa_ext_ch.*G_ampa_ext(pop_type);
                        s_nmda_ch = s_nmda(i-delay_ind, :);
                        s_gaba_ch = s_gaba(i-delay_ind, :);
                        if t(i) < t_task || t(i) > t_taskoff
                            I_temp = synapse_current(V_ch, s_ampa_ch, g_ampa_ext, s_nmda_ch, s_gaba_ch, ...
                            adja, AMPA, NMDA, GABA);
                        else
                            no_inhibition_adja = adja;
                            topN = 9;
                            no_inhibition_adja(N_E+1:end, 1:topN) = 0;
                            I_temp = synapse_current(V_ch, s_ampa_ch, g_ampa_ext, s_nmda_ch, s_gaba_ch, ...
                            no_inhibition_adja, AMPA, NMDA, GABA);
                        end
                        I_ch = sum(I_temp, 1);
                    end

                    %External input
                    spikes_ch = spikes(i, :);
                    if any(spikes_ch)
                        s_ampa_ext(i, spikes_ch) = s_ampa_ext(i, spikes_ch) + 1;
                    end

                    %Spiking Behavior
                    is_spike = (Vm(i, :)>= Vs);
                    if any(is_spike)
                        Vm(i, is_spike) = 0; %draw spike
                        Vm(i+1, is_spike) = Vr;
                        RP_ind(is_spike) = refract_ind(pop_type(is_spike));
                        s_ampa(i, is_spike) = s_ampa(i, is_spike) + 1;
                        x_nmda(i, is_spike) = x_nmda(i, is_spike) + 1;
                        s_gaba(i, is_spike) = s_gaba(i, is_spike) + 1;
                    end

                    %Synapse Update
                    s_ampa(i+1, :) = s_ampa(i, :) - dt*s_ampa(i, :) / tau_AMPA;
                    s_ampa_ext(i+1, :) = s_ampa_ext(i, :) - dt*s_ampa_ext(i, :)/tau_AMPA;
                    s_nmda(i+1, :) = s_nmda(i, :) + dt*(alpha*x_nmda(i, :).*(1-s_nmda(i, :)) ...
                        - s_nmda(i, :)/tau_NMDA_2);
                    x_nmda(i+1, :) = x_nmda(i, :) - dt*x_nmda(i, :)/tau_NMDA_1;
                    s_gaba(i+1, :) = s_gaba(i, :) - dt*s_gaba(i, :)/tau_GABA;

                    %Pulse Blocking Effects
                    if pulse
                        in_pp_rp = RP_pp_ind~=0;
                        in_ps_rp = RP_ps_ind~=0;
                        try
                            pulse_starters = I_ustim(i-1, :)==0 & I_ustim(i, :) ~= 0;
                        catch
                            assert(i==1, "Known indexing issue")
                            pulse_starters = false;
                        end
                        I_ustim_temp(i:i+stim_ind, in_pp_rp&pulse_starters) = 0;
                        I_ch(in_ps_rp) = 0;
                    end

                    %Voltage update (for non-refractory neurons)
                    non_rp = (RP_ind==0);
                    in_rp = (~non_rp)&(~is_spike);
                    dvdt = (-gL(pop_type).*(Vm(i, :)-EL) + I_ch + I_ustim_temp(i, :))./C(pop_type);
                    Vm(i+1, non_rp) = Vm(i, non_rp) + dvdt(non_rp)*dt;
                    Vm(i+1, in_rp) = Vr; %hold refractory neurons at Vr
                    RP_ind(in_rp) = RP_ind(in_rp) - 1; %decrement remaining refractory time

                    if pulse
                        %pulse refractory
                        try
                            pulse_enders = I_ustim(i-1, :) < 0 & I_ustim(i, :) == 0;
                            blocked_enders = pulse_enders & I_ustim_temp(i-1, :)==0;
                            unblocked_enders = pulse_enders & I_ustim_temp(i-1, :)~=0;
                            if any(blocked_enders)
                                RP_pp_ind(blocked_enders) = RP_pp_ind(blocked_enders) + ...
                                    get_blocked_RP(abs(I_ustim(i-1, blocked_enders)), I_b, t_pp, dt);
                                RP_ps_ind(blocked_enders) = RP_ps_ind(blocked_enders) + ...
                                    get_blocked_RP(abs(I_ustim(i-1, blocked_enders)), I_b, t_ps, dt);
                            end
                            if any(unblocked_enders)
                                RP_ps_ind(unblocked_enders) = get_RP(...
                                    abs(I_ustim(i-1, unblocked_enders)), I_b, t_ps, dt);
                                RP_pp_ind(unblocked_enders) = get_RP(...
                                    abs(I_ustim(i-1, unblocked_enders)), I_b, t_pp, dt);
                            end
                        catch
                            assert(i <= stim_ind, "Known Indexing Error")
                        end
                        RP_pp_ind(in_pp_rp) = RP_pp_ind(in_pp_rp) - 1;
                        RP_ps_ind(in_ps_rp) = RP_ps_ind(in_ps_rp) - 1;
                    end
                end 
                fast_parsave(output_trialpath, Vm);
            end
        end
    end
end
toc

function I = synapse_current(V, s_ampa, g_ampa_ext, s_nmda, s_gaba, adja, AMPA, NMDA, GABA)
    Mg = 1; %mM
    E_AMPA = 0e-3; %mV
    E_NMDA = 0e-3; %mV
    E_GABA = -70e-3; %mV
    
    %First, find the total synaptic conductance (summed over all synapses)
    %for each neuron (for AMPA, NMDA, and GABA)
    g_AMPA = s_ampa*(adja.*AMPA) + g_ampa_ext;
    g_NMDA = s_nmda*(adja.*NMDA);
    g_GABA = s_gaba*(adja.*GABA);
    
    %Then, use post-synaptic membrane potential to calculate total current
    I = zeros(3, length(V));
    I(1, :) = g_AMPA.*(E_AMPA-V);
    I(2, :) = g_NMDA.*(E_NMDA-V) ./ (1 + Mg*exp(-0.062*V*1000)/3.57);
    I(3, :) = g_GABA.*(E_GABA-V);
end

function RP_ind = get_RP(I_pulse, I_b, t_b, dt)
    RP = interp1(I_b, t_b, I_pulse);
    RP_ind = floor(RP / dt);
end

function RP_ind = get_blocked_RP(I_pulse, I_b, t_b, dt)
    RP = interp1(I_b, t_b, I_pulse) .* (I_pulse / I_b(end)).^2*2;
    RP_ind = floor(RP/dt);
end