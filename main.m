%%% Paul Adkisson
%%% 11.7.21
%%% Main Body of Wang (2002) Biophysical Attractor Model

clear;
sim_name = "P1_Rec";
sim_path = sprintf("Simulation %s", sim_name);
tic;
load(strcat(sim_path, "/bam_constants.mat"))
load(strcat(sim_path, "/adja.mat"))
load(strcat(sim_path, "/conductances.mat"))

for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j <= length(pulse_amps);
    if pulse
        fprintf("Pulse Stimulation Amplitude: %0.2fuA \n", stim_amp*1e6)
        input_stimpath = strcat(sim_path, sprintf("/ustim/%0.2fuA_pulse.mat", stim_amp*1e6));
        output_stimpath = strcat(sim_path, sprintf("/data/%0.2fuA_pulse", stim_amp*1e6));
        stim_coherences = pulse_coherences;
    elseif stim_amp ~= 0
        fprintf("Galvanic Stimulation Amplitude: %0.2fuA \n", stim_amp*1e6)
        input_stimpath = strcat(sim_path, sprintf("/ustim/%0.2fuA_galvanic.mat", stim_amp*1e6));
        output_stimpath = strcat(sim_path, sprintf("/data/%0.2fuA_galvanic", stim_amp*1e6));
        stim_coherences = galvanic_coherences;
    else
        disp("0nA Stimulation Control")
        input_stimpath = strcat(sim_path, sprintf("/ustim/%0.2fuA_galvanic.mat", stim_amp*1e6));
        output_stimpath = strcat(sim_path, sprintf("/data/%0.2fuA_galvanic", stim_amp*1e6));
        stim_coherences = control_coherences;
    end
    load(input_stimpath, 'I_ustim', 'Vmir')
    mkdir(output_stimpath)
    for c = stim_coherences
        fprintf("Coherence: %0.1f%% \n", c*100)
        input_coherentpath = sprintf("Simulation %s/spikes/c=%0.3f", sim_name, c);
        output_coherentpath = strcat(output_stimpath, sprintf("/c=%0.3f", c));
        mkdir(output_coherentpath)
        parfor trial = start_trial:end_trial
        %for trial = start_trial:end_trial
            fprintf("trial: %0.0f \n", trial)
            input_trialpath = strcat(input_coherentpath, sprintf("/trial%0.0f.mat", trial));
            output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));

            %external input
            input_spikes = load(input_trialpath, 'spikes');
            spikes = input_spikes.spikes;

            %ustim
            I_ustim_temp = I_ustim;
            pblocked = zeros(1, N); %neurons blocked at the start of each pulse

            %LIF variables
            RP_ind = zeros(1, N);
            if pulse
                RP_pp_ind = zeros(1, N);
                RP_ps_ind = zeros(1, N);
            end
            Vm = EL*ones(length(t), N);
            I_AMPA_ext = zeros(length(t), 1);
            I_AMPA_rec = zeros(length(t), 1);
            I_NMDA = zeros(length(t), 1);
            I_GABA = zeros(length(t), 1);
            I_LEAK = zeros(length(t), 1);
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
                    I_temp = synapse_current(V_ch, s_ampa_ch, g_ampa_ext, s_nmda_ch, s_gaba_ch, ...
                        adja, AMPA, NMDA, GABA);
                    I_ch = sum(I_temp, 1);
                    % Record currents
                    I_AMPA_ext(i) = mean(g_ampa_ext(1:num_group).*(0-Vm(i, 1:num_group)));
                    I_AMPA = mean(I_temp(1, 1:num_group));
                    I_AMPA_rec(i) = I_AMPA - I_AMPA_ext(i);
                    I_NMDA(i) = mean(I_temp(2, 1:num_group));
                    I_GABA(i) = mean(I_temp(3, 1:num_group));
                    I_LEAK(i) = mean(-gL(1).*(Vm(i, 1:num_group)-EL));
                end

                %External input
                spikes_ch = spikes(i, :);
                if any(spikes_ch)
                    s_ampa_ext(i, spikes_ch) = s_ampa_ext(i, spikes_ch) + 1;
                end

                %Spiking Behavior
                if pulse
                    in_pp_rp = RP_pp_ind~=0;
                    in_ps_rp = RP_ps_ind~=0;
                    is_spike = (Vm(i, :)>= Vs) & (~in_ps_rp | (I_ustim(i, :)>0 & ~in_pp_rp));
                else
                    is_spike = (Vm(i, :)>= Vs);
                end
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

                if pulse
                    %Pulse Blocking Effects
                    try
                        pulse_starters = I_ustim(i-1, :)==0 & I_ustim(i, :) ~= 0;
                        pulse_flippers = I_ustim(i-1, :)~=0 & I_ustim(i-1, :)==-I_ustim(i, :);
                        pulse_enders = I_ustim(i-1, :) < 0 & I_ustim(i, :) == 0;
                    catch
                        assert(i==1, "Known indexing issue")
                        pulse_starters = false;
                        pulse_flippers = false;
                        pulse_enders = false;
                    end
                    I_ustim_temp(i:i+stim_ind, in_pp_rp&pulse_starters) = 0;
                else
                    %Depolarization Block
                    I_tot = -gL(pop_type).*(Vm(i, :)-EL) + I_ch + I_ustim_temp(i, :);
                    depol_blocked = I_tot > depol_block_thresh;
                    I_ch(depol_blocked) = 0;
                    I_ustim_temp(i, depol_blocked) = 0;
                end

                %Voltage update (for non-refractory neurons)
                non_rp = (RP_ind==0);
                in_rp = (~non_rp)&(~is_spike);
                if pulse
                    if any(pulse_starters)
                        pblocked(pulse_starters) = in_rp(pulse_starters) | in_pp_rp(pulse_starters);
                    end
                    pstart = pulse_starters&~pblocked;
                    pflip = pulse_flippers&~pblocked&non_rp;
                    pend = pulse_enders&~pblocked&non_rp;
                    %mirror pulse estimate
                    Vm(i, pstart) = Vm(i, pstart) + Vmir(pstart);
                    Vm(i, pflip) = Vm(i, pflip) - 2*Vmir(pflip);
                    Vm(i, pend) = Vm(i, pend) + Vmir(pend);
                end
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
                            %neurons with pulse amplitudes large enough to experience residual pulse-pulse blocking
                            pp_blockers = abs(I_ustim(i-1, :)) > I_b(min_pp_idx) & blocked_enders;
                            RP_pp_ind(pp_blockers) = get_RP(abs(I_ustim(i-1, pp_blockers)), I_b, t_pp, dt) / 2;
                            RP_ps_ind(blocked_enders) = get_RP(abs(I_ustim(i-1, blocked_enders)), I_b, t_ps, dt);
                        end
                        if any(unblocked_enders)
                            RP_ps_ind(unblocked_enders) = get_RP(...
                                abs(I_ustim(i-1, unblocked_enders)), I_b, t_ps, dt);
                            RP_pp_ind(unblocked_enders) = get_RP(...
                                abs(I_ustim(i-1, unblocked_enders)), I_b, t_pp, dt);
                        end
                    catch
                        assert(i <= stim_ind, "Pulse Refractory Error")
                    end
                    %RP's cannot go below 0
                    RP_pp_ind(in_pp_rp) = max(RP_pp_ind(in_pp_rp) - 1, 0);
                    RP_ps_ind(in_ps_rp) = max(RP_ps_ind(in_ps_rp) - 1, 0);
                end
            end 
            fast_parsave(output_trialpath, Vm);
            % save(output_trialpath, 'Vm', 'I_AMPA_ext', 'I_AMPA_rec', ...
            %     'I_NMDA', 'I_GABA', 'I_LEAK')
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